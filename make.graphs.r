rm(list=ls())
graphics.off()

####################################################################################################
### FUNCTION DEFINITIONS
####################################################################################################

fix.subject.ids <- function (in.subject.ids) {

    in.subject.ids=gsub("_[ABCD]", "", as.character(in.subject.ids), fixed=FALSE)
    subjectIds=as.character(in.subject.ids)
    if (length(pmatch(subjectIds, "300")) > 0 ) {
        subjectIds=sub("^300", "169/300", subjectIds, fixed=FALSE)
        subjectIds=sub("^0", "", subjectIds, fixed=FALSE)
    }
    return (subjectIds)
}

drop.subjects <- function (in.netcc.filenames.and.subjects, in.drop.subject.list) {
    if (length(in.drop.subject.list) > 0 ) {
        df=in.netcc.filenames.and.subjects[ ! in.netcc.filenames.and.subjects$Study.ID %in% in.drop.subject.list, ]
    }
    return(df)
}

read.demographics.table <- function(in.filename) {
    cat("*** Reading", in.filename, "\n")
    demographics=read.csv(in.filename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
    cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

    return(demographics)
}

fix.demographics.table <- function(in.demographics) {
    in.demographics=in.demographics[grep("Test_Subject", in.demographics$ID, invert=TRUE), ]
    
    ## in the case where the ID does not contain a _ followed by a
    ## timepoint letter separate will emit a warning to indicate that
    ## it got two few values for the second component (the timepoint
    ## below). Use supressWarnings to suppress these warnings
    in.demographics = suppressWarnings(separate_(in.demographics, "ID", into=c("Study.ID", "timepoint"), sep="_", remove=TRUE))
    if (any(grepl("169/300", in.demographics$Study.ID)) ) {
        cat("*** Replacing subject ID 169/300 with 300\n")        
        in.demographics$Study.ID=as.factor(sub("^169/300", "300", as.character(in.demographics$Study.ID), fixed=FALSE))
    }

    if (any(grepl("378", in.demographics$Study.ID)) ) {
        cat("*** Setting gender for subject 378\n")
        ## 378 identified as transgender, was biologically female and not on any hormone therapy
        in.demographics[in.demographics$Study.ID=="378", "Gender"]="F"
        in.demographics$Gender=droplevels(in.demographics$Gender)
    }

    cat("*** Renaming Grp to Group\n")
    in.demographics = rename(in.demographics, replace=c("Grp" = "Group"))
    
    return(in.demographics)
}

read.wasi.table <- function(in.filename) {
    cat("*** Reading", in.filename, "\n")
    wasi=read.csv(in.filename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
    cat(sprintf("*** Read WASI-II data for %s unique subjects\n",  length(unique(wasi$SubID))))

    return(wasi)
}

fix.table <- function (in.data, in.id.column.name=NULL) {

    if (is.null(in.id.column.name)) {
        stop("*** in.id.column.name cannot be NULL. it must specify a vaild column name\n")
    }

    ## in the case where the ID does not contain a _ followed by a
    ## timepoint letter separate will emit a warning to indicate that
    ## it got two few values for the second component (the timepoint
    ## below). Use supressWarnings to suppress these warnings
    in.data = suppressWarnings(separate_(in.data, in.id.column.name, into=c("Study.ID", "timepoint"), sep="_", remove=TRUE))
    
    if (any(grepl("169/300", in.data$Study.ID)) ) {
        cat("*** Replacing subject ID 169/300 with 300\n")
        in.data$Study.ID=as.factor(sub("^169/300", "300", as.character(in.data$Study.ID), fixed=FALSE))
    }

    if (any(grepl("169", in.data$Study.ID)) ) {
        cat("*** Replacing subject ID 169 with 300\n")
        in.data$Study.ID=as.factor(sub("^169", "300", as.character(in.data$Study.ID), fixed=FALSE))
    }

    return(in.data)
}

read.telomere.table <- function(in.filename) {
    cat("*** Reading", in.filename, "\n")
    bio.data=read.csv(in.filename, na.strings=c("degraded", "T/S too high", "not enough DNA"), header=TRUE)

    return(bio.data)
}

read.mt.dna.table = read.telomere.table


fix.dates <- function (inData) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$MRI)
    
    ## now convert to year/month/day
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$MRI)
    
    inData$DOB=as.Date(inData$DOB, "%y/%m/%d")
    inData$MRI=as.Date(inData$MRI, "%y/%m/%d")

    return(inData)
}

compute.age <- function(inData) {

    age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inData$age.in.years=age.in.weeks/52

    return(inData)
}

build.netcc.filenames <- function (in.subjects, in.graph.subdir) {

    filenames=sapply(in.subjects,
        function(x) {
            sprintf("%s/%s/%s/%s.pm.cleanEPI.aal2.whole.ts.01_000.netcc", data.dir, x, in.graph.subdir, x, x)
        })
    
    return(filenames)
}


check.for.motion.exclusion.files <- function (in.subjects) {

    build.do.not.analyze.filenames <- function(in.subjects, in.threshold=20) {
        filenames=sapply(in.subjects,
            function(ss) {
                sprintf("%s/%s/rsfcPreprocessed/00_DO_NOT_ANALYSE_%s_%dpercent.txt",  data.dir, ss, ss, in.threshold)
            })
        
        return(filenames)
    }
    
    df = build.do.not.analyze.filenames(in.subjects)
    df = data.frame("do.not.analyze.filename"=df, "excessive.motion"=file.exists(df))
    df$ID=rownames(df)

    drop.subject.list=c()
    motion.contaminated.subjects.count=sum(df$excessive.motion)
    if (motion.contaminated.subjects.count > 0) {

        drop.subject.list=df[df$excessive.motion==TRUE, "ID"]

        cat("***", motion.contaminated.subjects.count, "of", length(in.subjects),
            paste("(", round(motion.contaminated.subjects.count/length(in.subjects) * 100 , 2) ,"%)", sep=""),  "are excessively contaminated by motion\n")
        cat("*** The following subjects are excessively contaminated by motion\n")
        cat("---", str_wrap(paste (drop.subject.list, collapse=" "), width=80), "\n")

        cat("*** The following subjects are NOT excessively contaminated by motion\n")
        cat("+++", str_wrap(paste (df[df$excessive.motion==FALSE, "ID"], collapse=" "), width=80), "\n")

    }
    
    return(drop.subject.list)
}

filter.nonexistant.files <- function (in.file.list, in.print.nonexistant=FALSE) {
    existant=sapply(in.file.list, file.exists)
    if (in.print.nonexistant && sum(!existant) > 0) {
        cat("*** The following netcc files do not exist:\n")
        cat(paste("---", in.file.list[!existant]), sep="\n")
    }
    
    aa=in.file.list[existant]
    df = data.frame("filename"=aa, "subject.name"=names(aa))
    df = df %>% separate(subject.name, into=c("Study.ID", "timepoint"), sep="_", remove=FALSE)

    return(df)
}

read.netcc.files <- function (in.netcc.filenames, in.atlas.dt, in.matrix.type="CC") {
    dnames=list()
    dnames[[1]] = in.atlas.dt$name
    dnames[[2]] = in.atlas.dt$name
    dnames[[3]] = in.netcc.filenames$Study.ID
    ## print(dnames)
    
    netcc=array(as.numeric(NA),
        dim=c(length(dnames[[1]])[1], length(dnames[[2]])[1], length(in.netcc.filenames$subject.name)),
        dimnames=dnames)
    ## cat("*** netcc array is", paste(dim(netcc)), "\n")
    for (ii in seq.int(1, length(in.netcc.filenames$subject.name))) {

        lines=readLines(as.character(in.netcc.filenames[ii, "filename"]))
        matrix.header.line.number=grep(paste("^# ", in.matrix.type, "$", sep=""), lines)
        if (length(matrix.header.line.number) == 0 ) {
            stop(paste("Cannot find matrix type", in.matrix.type, "in", in.netcc.filenames[ii, "filename"], "\n"))
        }

        cat(sprintf("+++ Reading (%04d/%04d) (matrix %s begins at line %04d) %s\r", ii, length(in.netcc.filenames$subject.name),
                    in.matrix.type, matrix.header.line.number, in.netcc.filenames[ii, "filename"]))

        ff=read.table(textConnection(lines), header=FALSE, skip=matrix.header.line.number, nrows=length(dnames[[1]]))        
        ## print(dim(ff))
        netcc[, , ii ] = as.matrix(ff)
    }
    cat("\n")
    return(netcc)
}

check.for.zero.matrix.subjects <- function(in.netcc) {

    zero.subjects=vector(mode="numeric", length=dim(in.netcc)[3])
    names(zero.subjects)=dimnames(in.netcc)[[3]]

    for (ss in seq.int(from=1, to=dim(in.netcc)[3], by=1)) {
        zero.subjects[ss] = all(in.netcc[ , , ss] == 0)
    }

    drop.subject.list=c()
    zero.count.subjects=sum(zero.subjects)
    if (zero.count.subjects > 0) {

        drop.subject.list=names( zero.subjects)[ zero.subjects == TRUE ]
        cat("***", zero.count.subjects, "of", dim(in.netcc)[3], paste("(", round(zero.count.subjects/dim(in.netcc)[3] * 100 , 2) ,"%)", sep=""),  "have all zero correlation matrices\n")
        cat("*** The following subjects have all zero correlation matrices\n")
        cat("---", str_wrap(paste (drop.subject.list, collapse=" "), width=80), "\n")

        cat("*** The following subjects have NON zero correlation matrices\n")
        cat(str_wrap(paste (names( zero.subjects)[ zero.subjects == FALSE ], collapse=" "), indent=4, exdent=4, width=80), "\n")

    }

    return(drop.subject.list)
}


drop.subjects.from.netcc <- function (in.netcc, in.drop.subject.list){
    if(length(in.drop.subject.list) > 0) {
        ## cat("*** Dropping subjects with all zero correlation matrices\n")
        return(in.netcc[, , -(which (names(in.netcc[1, 1, ]) %in% in.drop.subject.list))])
    } else {
        return(in.netcc)
    }
}

convert.to.feature.matrix <- function(in.netcc, in.correlation.tag) {

    ut=upper.tri(in.netcc[, , 1])
    ## upper.tri privides a matrix of TRUE/FALSE values. When this is
    ## summed (as below) it gives the total number of elements to be
    ## retained from each slice of in.netcc and hence the number of
    ## columns (features) in fv per subject
    fv=matrix(in.netcc[ut],
        nrow=dim(in.netcc)[3], ## number of subjects in the netcc 3D array
        ncol=sum(ut), byrow=TRUE)

    ## this gives us the array indices of the elements in the upper
    ## triangle. 
    nn=which(ut, arr=TRUE)
    ## We can then use these pairs of indices to concatentate the
    ## corresponding elements from the rownames and colnames of
    ## in.netcc to provide feature vector names that preserve the
    ## pairwise nature of the correlation coefficients stored in the
    ## in.netcc matrix thus allowing us to know to what pair of brain
    ## regions each element in the feature vector corresponds
    rn=rownames(in.netcc)
    cn=colnames(in.netcc)
    ## the labels will be of the form R000.R000 which can be
    ## translated back to ROI names by using the integer component to
    ## index the aal labels data frame
    colnames(fv)=apply(nn, 1, function(xx) { sprintf("%s.%s.%s", in.correlation.tag, rn[xx[1]], cn[xx[2]]) } )
    rownames(fv)=names(in.netcc[1, 1, ])

    return(fv)
}

list.excluded.subjects <- function(in.df) {

    orig.dims=dim(in.df)
    dl=c()
    print.dropped.subject.list <- function (in.drop.list) {
        if (any(in.drop.list)) {
            cat ("--- The following subjects are to dropped from the data frame\n")
            dl=as.vector (  in.df[which(in.drop.list), "Study.ID"])
            cat ("---", paste (dl, collapse=" "), "\n")
            cat ("---", paste (as.vector (  in.df[which(in.drop.list), "Group"]),     collapse=" "), "\n")
        }
        return(dl)
    }
    
    drop.count=0
    ## cat("*** Filtering to exlcude subjects who were medicated or do not meet other inclusion/exclusion criteria\n")
    ## cat("*** Subject numbers BEFORE filtering\n")
    ## print(addmargins(table(in.df[, c("Gender", "Group")])))

    ## ################################################
    ## Drop MDDs with CDRS-R > 54
    cat("*** Checking for NCLs with CDRS-R tscore > 54\n")
    drop.list.cdrsr.gt.fiftyfour=in.df$Group=="NCL" & in.df$CDRS.tscore > 54
    drop.count=drop.count + sum(drop.list.cdrsr.gt.fiftyfour)
    dropped.subjects.list=print.dropped.subject.list(drop.list.cdrsr.gt.fiftyfour)
    ## in.df = in.df[ ! drop.list , ]

    ## ################################################    
    ## Drop NCLs with CDRS-R > 55
    cat("*** Checking for MDDs with CDRS-R tscore < 55\n")
    drop.list.cdrsr.lt.fiftyfive=in.df$Group=="MDD" & in.df$CDRS.tscore < 55
    drop.count=drop.count + sum(drop.list.cdrsr.lt.fiftyfive)    
    dropped.subjects.list=c(dropped.subjects.list, print.dropped.subject.list(drop.list.cdrsr.lt.fiftyfive))
    ## in.df = in.df[ ! drop.list, ]    

    ## ################################################    
    ## Drop medicated subjects (This list was compiled for the TL
    ## paper by Laura. See
    ## rsfcGraphAnalysis/scripts/analyseGroupSubjects.r and
    ## structuralTelomereAnalysis/scripts/analyseStructuralVolumes_Laura.r
    ## Now remove subjects that are NOT medication-naive (Prozac,
    ## Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbalta and
    ## ADHD meds)
    cat("*** Checking for subjects known to have been on medications\n")
    drop.list.medicated=in.df$Study.ID %in%
    c("130", "132", "318", "319", "320",
      "322", "323", "324", "325", "329", "345", "376", "384", "333", "349")
    drop.count=drop.count + sum(drop.list.medicated)
    dropped.subjects.list=c(dropped.subjects.list, print.dropped.subject.list(drop.list.medicated))
    ## in.df = in.df[ ! drop.list, ]

    ## ################################################    
    ## Drop subjects with WASI Performance and WASI Verbal < 80
    cat("*** Checking for WASI Verbal < 80 and WASI Perfromance < 80\n")
    drop.list.wasi=in.df$Verbal < 80 & in.df$Performance < 80
    drop.count=drop.count + sum(drop.list.wasi)
    dropped.subjects.list=c(dropped.subjects.list, print.dropped.subject.list(drop.list.wasi))
    ## in.df = in.df[ ! drop.list, ]

    ## 391 is in group RUM on the demographics sheet, but I've no idea
    ## what this means so dump them
    cat("*** Checking for other excluded subjects\n")
    drop.list.rum=in.df$Study.ID %in% c("391")
    drop.count=drop.count + sum(drop.list.rum)
    dropped.subjects.list=c(dropped.subjects.list, print.dropped.subject.list(drop.list.rum))
    ## in.df = in.df[ ! drop.list, ]

    cat("***", drop.count, "subjects are to be dropped\n")
    ## cat("*** Subject numbers AFTER filtering\n")
    ## print(addmargins(table(in.df[, c("Gender", "Group")])))

    ## new.dims=dim(in.df)
    ## if (isTRUE(all.equal(orig.dims, new.dims))) {
    ##     cat("*** No subjects were filtered out\n")
    ## }
    
    return(dropped.subjects.list)
}


create.heatmap <- function(in.netcc) {

    library(reshape2)
    library(RColorBrewer)
    library(ggplot2)
    library(ggdendro)

    ## netcc.mean=apply(netcc[1:30, 1:30, ], c(1, 2), mean, na.rm=TRUE)
    netcc.mean=apply(netcc, c(1, 2), mean, na.rm=TRUE)    

    ## brewer.RdBu=rev(brewer.pal(11, "RdBu"))
    ## color.palette = colorRampPalette(brewer.RdBu, space = "Lab")

    dd.col=as.dendrogram(hclust(dist(netcc.mean)))
    col.ord=order.dendrogram(dd.col)

    ## Use correlation between variables as distance
    ## dd <- as.dist((1-netcc.mean)/2)
    ## dd=dist(netcc.mean)
    dd.row=as.dendrogram(hclust(dist(netcc.mean)))
    row.ord=order.dendrogram(dd.row)
    
    xx=netcc.mean[col.ord, row.ord]
    xx_names = attr(xx, "dimnames")
    df = as.data.frame(xx)
    colnames(df) = xx_names[[2]]
    df$roi=xx_names[[1]]
    df$roi=with(df, factor(roi, levels=roi, ordered=TRUE))
    
    mdf = melt(df, id.vars="roi")
    ddata_x <- dendro_data(dd.row)
    ddata_y <- dendro_data(dd.col)
    
    my.base.size=8
    my.heat.theme=
        theme_bw(base_size =  my.base.size) +
            theme(
                ## legend.position="none",
                legend.position="bottom",        
                ## panel.grid.major = element_blank(),
                ## panel.grid.minor = element_blank(),
                
                ##remove the panel border
                ## panel.border = element_blank(),
                
                ## add back the axis lines
                axis.line=element_line(colour = "grey50"),
                
                ##axis.title.x=element_blank(),
                axis.title.x = element_text(size=my.base.size, vjust=0),
                axis.text.x = element_text(size=my.base.size, angle=90),
                axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
                plot.title=element_text(size=my.base.size*1.2, vjust=1))
    
    heatmap <- ggplot(mdf, aes(x=variable, y=roi))
    heatmap=heatmap + geom_tile(aes(fill=value))
    ## heatmap=heatmap + scale_fill_gradient2(name="Correlation", low="#053061", mid="#F7F7F7", high="#67001F", limit = c(min(min(netcc.mean)), max(max(netcc.mean))))
    heatmap=heatmap + scale_fill_gradient2(name="Correlation", low="#053061", mid="#F7F7F7", high="#67001F", limit = c(-1, 1))    
    heatmap=heatmap + labs(x="ROI", y="ROI")
    heatmap=heatmap + my.heat.theme
    

    return(heatmap)
}

filter.near.zero.variance.predictors <- function (in.df, in.predictor.tag=NULL) {
    if (! require(caret)) {
        stop("*** Couldn't load the caret library\n")
    }

    if (is.null(in.predictor.tag)) {
        stop("*** No predictor variable tag specified. Stopping\n")
    }

    nzv = nearZeroVar(feature.df[, grep(in.predictor.tag, colnames(feature.df), fixed=TRUE)], saveMetrics= TRUE)

    if (any(nzv$nzv)) {
        cat("*** The following predictors have near zero variance and will be removed:\n")
        cat ("---", paste (rownames(nzv)[nzv$nzv], collapse=" "), "\n")
    }

    ## return all columns but those with near zero variance
    return(in.df[,  ! (colnames(in.df) %in% rownames(nzv)[nzv$nzv])])
}

check.number.of.rois.match <- function () {
    ## check that the number of ROIs in the list of ROI names matches the
    ## number of ROI on the rows and columns dimensions of the netcc array
    if ( ! isTRUE(all.equal(rep( dim(atlas.dt)[1], 2), dim(netcc)[1:2])) ) {
        stop("*** The number of ROIs in the atlas (", dim(atlas.dt)[1],
             ") does not match the number of ROIs on in the first two dimensions (",
             paste(dim(netcc)[1:2], collapse=", "),
             ") of the netcc array\n", sep="")
    }
}


## zscore.feature.matrix <- function(in.feature.matrix) {
##     cat("*** Z-scoring feature matrix\n")
##     ## atanh is equivalent to Fisher's r-to-z transform
##     apply(in.feature.matrix, 2, atanh)
## }

filter.feature.matrix <- function(in.feature.matrix, in.p.threshold=0.01) {
    cat("*** Filtering feature matrix\n")
    less.than.threshold=apply(feature.matrix, 2, function (xx) { t.test(xx)$p.value < in.p.threshold } )
    count.below.threshold=sum(less.than.threshold)
    if(count.below.threshold > 0 ) {
        cat("***", count.below.threshold, "of", dim(in.feature.matrix)[2],
            paste("(", round(count.below.threshold/dim(in.feature.matrix)[2] * 100 , 2) ,"%)", sep=""),  "are below the", in.p.threshold, "p threshold on the t tests against 0\n")
        in.feature.matrix = in.feature.matrix[, less.than.threshold]
    }
    return(in.feature.matrix)
}


threshold.netcc <- function (in.netcc.matrix, in.density=0.1) {
    N <- ncol(in.netcc.matrix)
    emax <- N  * (N - 1) / 2

    out=list()

    ## the following code is from corr_matrix.R in brainGraph
    ##
    ## try to make the data structures as similar as possible, if
    ## not identical, so that we can use brainGraph for other
    ## processess later if necessary
    thresh=sort(in.netcc.matrix[lower.tri(in.netcc.matrix)])[emax - in.density * emax]
    r.thresh <- ifelse(in.netcc.matrix > thresh, 1, 0)
    weights=in.netcc.matrix * r.thresh
    diag(weights) = 0
    
    out <- list(R=in.netcc.matrix, P=pnorm(in.netcc.matrix), r.thresh=r.thresh, weights=weights, threshold=thresh)
    
    return(out)
}


threshold.correlations <- function (in.netcc, in.graph.densities) {
    tc =
        foreach(sn=icount(dim(in.netcc)[3])) %dopar% {
            foreach(gd=in.graph.densities) %do% {
                threshold.netcc(in.netcc[, , sn], gd)
            }}
    ## now set the names
    tc=lapply(tc, function (xx) { names(xx) = in.graph.densities; return(xx) })
    names(tc)=names(in.netcc[1, 1, ])
    
    return(tc)
}

## ##################################################################################################
## END OF FUNCITONS
## ##################################################################################################

### setup path variables
if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    max.cpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    max.cpus=as.integer(system("getconf _NPROCESSORS_ONLN", intern=TRUE))
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
    max.cpus=NA
}

source("common.functions.r")

library(stringr)
library(tidyr)
library(plyr)
library(data.table)
library(foreach)
library(iterators)

study.root.dir=file.path(root.dir, "sanDiego/machLearnT1Analysis")
standard.data.dir=file.path(study.root.dir, "standard")

scripts.dir=file.path(study.root.dir, "scripts")
data.dir=file.path(study.root.dir, "data")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
group.results.dir=file.path(data.dir, "Group.results")

subjects=dir(data.dir, pattern="[0-9][0-9][0-9]_A*")

## before we go any further we need to deal with the issue of subject 169/300 having two IDs
## this subject is known as 169 in wasi.data, telomere.data, and mt.dna.data
## 169/300 in the demographics file, and 300 is used for their MRI data.
##
## To simplify things we'll standardize on using 300 since it will
## match the namd used for the MRI data which is considerably less
## easy to change

## we also take the opportunity to standardize on the use of Study.ID
## as the name of the column that identidye subjects

demographics.filename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=fix.demographics.table(read.demographics.table(demographics.filename))

telomere.filename=file.path(admin.data.dir, "Tony Yang TL data 06162014.csv")
telomere.data=fix.table(read.telomere.table(telomere.filename), "sample.ID")

mt.dna.filename=file.path(admin.data.dir, "mtDNA.csv")
mt.dna.data=fix.table(read.mt.dna.table(mt.dna.filename), "subject") 

wasi.filename=file.path(admin.data.dir, "WASI.csv")
wasi.data=fix.table(read.wasi.table(wasi.filename), "SubID")

####################################################################################################
## Now build the lists of netcc filenames and filter out subjects
## contaminated by motion, that do not meet inclusion/exclusion
## criteria, have all zero correlation matrices, 

## now read the graphs created from the RSFC analysis (note the
## "rsfcGraphs" argument) netcc filenames that do not exist are
## filtered from the list


library(devtools)
## library(brainGraph)
devtools::load_all(file.path(Sys.getenv("HOME"), "src", "brainGraph"))

atlas="brainGraph::aal2.94"
atlas.dt=eval(parse(text=atlas))
cat("*** Using the", atlas, "atlas\n")

## use the partial correlation (PC) coefficients created by 3dNetCorr
matrix.type="PC"
force.netcc.generation=FALSE
saved.netcc.filename=file.path(group.results.dir, paste("netcc", matrix.type, "Rdata", sep="."))
if (file.exists(saved.netcc.filename) && ! force.netcc.generation) {
    cat("*** Loading pregenerated netcc array and netcc.filenames.and.subjects.df from", saved.netcc.filename, "\n")
    load(saved.netcc.filename)

    cat(sprintf("*** Loaded saved netcc array with %d rows, %d columns, and %d slices (subjects)\n",
                dim(netcc)[1], dim(netcc)[2], dim(netcc)[3]))

} else {

    netcc.filenames.and.subjects.df=filter.nonexistant.files(build.netcc.filenames(subjects, "rsfcGraphs"), in.print.nonexistant=TRUE)

    ## check for subjects that should be excluded because the moved too
    ## much during the scan
    ## motion.exclusion.subject.list=check.for.motion.exclusion.files(subjects)
    ## netcc.filenames.and.subjects=drop.subjects(netcc.filenames.and.subjects,
    ## motion.exclusion.subject.list)
    
    ## build the lists of netcc filenames and read them in
    netcc=read.netcc.files(netcc.filenames.and.subjects.df, atlas.dt, in.matrix.type=matrix.type)
    cat("*** Saving netcc array and netcc.filenames.and.subjects.df to", saved.netcc.filename, "\n")
    pigz.save(netcc, netcc.filenames.and.subjects.df, file=saved.netcc.filename)
}
## check for subjects with all zero correlation matrices and then drop
## them
drop.subject.list=check.for.zero.matrix.subjects(netcc)
## print(drop.subject.list)
## print(dim(netcc))
netcc=drop.subjects.from.netcc(netcc, drop.subject.list)
## print(dim(netcc))
## print(dim(netcc.filenames.and.subjects.df))
netcc.filenames.and.subjects.df=drop.subjects(netcc.filenames.and.subjects.df, drop.subject.list)
## print(dim(netcc.filenames.and.subjects.df))

## check that the number of ROIs in the netcc files (i.e., the number
## of columns and rows since the matrices should be square), match the
## number of ROIs read in from the table that contains thre names of
## the ROIs in the AAL atlas used to generate the netcc files in the
## first place
check.number.of.rois.match()

## hm=create.heatmap(netcc)
## print(hm)


####################################################################################################
## feature.matrix=convert.to.feature.matrix(netcc, "RSFC")
## feature.matrix=zscore.feature.matrix(feature.matrix)
## feature.matrix=filter.feature.matrix(feature.matrix, in.p.threshold=0.01)
## feature.matrix=convert.to.feature.matrix(netcc[1:30, 1:30, 1:30], "RSFC")

saved.graph.data.structures.filename=file.path(group.results.dir, paste("graphs", matrix.type, "Rdata", sep="."))
## a list of data structures to be saved so that can be reloaded
## later to obviate the need to renenerate all the graphs and
## associated metrics
save.structures.list=list()

characteristics.df=data.frame(
    "Study.ID"           =netcc.filenames.and.subjects.df$Study.ID,
    demographics                      [match(netcc.filenames.and.subjects.df$Study.ID, demographics$Study.ID),  c("Group", "Gender", "DOB", "MRI", "CDRS.tscore")],
    wasi.data                         [match(netcc.filenames.and.subjects.df$Study.ID, wasi.data$Study.ID),     c("Verbal", "Performance", "Full")],
    "telomere.final.T.S"=telomere.data[match(netcc.filenames.and.subjects.df$Study.ID, telomere.data$Study.ID), "final.T.S"],
    "mt.dna"            =mt.dna.data  [match(netcc.filenames.and.subjects.df$Study.ID, mt.dna.data$Study.ID),   "mtDNA"])

characteristics.df=fix.dates(characteristics.df)
characteristics.df=compute.age(characteristics.df)

## drop these two columns as they are no longer needed
characteristics.df$DOB=NULL
characteristics.df$MRI=NULL
og.characteristics.df = characteristics.df

## list medicated subjects, subjects with too high/low CDRS-R and WASI score
drop.subject.list=list.excluded.subjects(characteristics.df)
## print(drop.subject.list)
## drop the subjects from the characteristics df
## cat("*** BEFORE removing subjects from characteristics.df", dim(characteristics.df), '\n')
characteristics.df=characteristics.df[ ! characteristics.df$Study.ID %in% drop.subject.list, ]
## cat("*** AFTER removing subjects from characteristics.df", dim(characteristics.df), '\n')
characteristics.df=droplevels(characteristics.df)

save.structures.list=append(save.structures.list, "characteristics.df")

## drop.subject.list=paste(drop.subject.list, "_A", sep="")
## remove subjects from the netcc array
## cat("*** BEFORE removing subjects from netcc", dim(netcc), '\n')
netcc=drop.subjects.from.netcc(netcc, drop.subject.list)
## print(dim(netcc))
## cat("*** AFTER removing subjects from netcc", dim(netcc), '\n')

## remove subjects from the netcc.filenames.and.subjects.df
## cat("*** BEFORE removing subjects from netcc.filenames.and.subjects.df", dim(netcc.filenames.and.subjects.df), '\n')
netcc.filenames.and.subjects.df=drop.subjects(netcc.filenames.and.subjects.df, drop.subject.list)
save.structures.list=append(save.structures.list, "netcc.filenames.and.subjects.df")
## cat("*** AFTER removing subjects from netcc.filenames.and.subjects.df", dim(netcc.filenames.and.subjects.df), '\n')

## a final sanity check to ensure the subject IDs all match in all of
## the data frames and arrays
aa=cbind(as.character(characteristics.df$Study.ID), names(netcc[1, 1, ]), as.character(netcc.filenames.and.subjects.df$Study.ID))
if (! all(apply(aa, 1, function (xx) { ! length(unique(xx)) > 1 } ))) {
    stop("*** The subjects IDs in charactersitics.df, netcc[1, 1, ], and netcc.filenames.and.subjects.df$Study.ID do not match.\nCannot continue\n")
} else {
    cat("*** Subject names in charactersitics.df, netcc, and netcc.filenames.and.subjects.df all match\n")
}
## delete aa, it's no longer needed
rm(aa)

## ##############################################################################
## now that all of the filtering has been done it's time to z-score the array
cat("*** Applying Fisher's r-to-z scoring to the netcc array\n")
netcc.og=netcc
netcc=atanh(netcc)
save.structures.list=append(save.structures.list, "netcc")

## ##############################################################################
cat("*** The subject distribution is as follows:\n")
subject.distribution=addmargins(xtabs(~ Group + Gender, data=characteristics.df))
print(subject.distribution)

## ##############################################################################
## if you only want to use one density set that here to a vector of
## length one
graph.densities=densities=seq.int(0.1, 0.6, 0.01)
## graph.densities=densities=c(0.05)
save.structures.list=append(save.structures.list, "graph.densities")
cat("*** Creating graphs for each subject at the following",
    length(graph.densities),
    ifelse(length(graph.densities) > 1, "densities:", "density:"), "\n",
    str_wrap(paste (graph.densities, collapse=" "), indent=3, exdent=4, width=80) ,"\n")

## ##############################################################################
parallel.executation=TRUE
if (parallel.executation) {
    cat("*** Enabling parallel processing\n")
    library(doMC)
    ## why 51? because there are 51 densities set up below. It's only
    ## choosen to get as many densities processed in parallel as
    ## possible
    registerDoMC(cores = min(51, max.cpus))
    cat("*** Using", getDoParWorkers(), "CPU cores\n")
    cat("*** Plyr progress bars are disabled when parallel computation is enabled\n")
    progress.bar.type='none'        
} else {
    cat("*** Plyr progress bars are set to text\n")        
    progress.bar.type='text'
}

## ##############################################################################
cat("*** Thresholding correlations at", length(graph.densities),
    ifelse(length(graph.densities) > 1,
           "different densities", "density"), "for each subject\n")

## this commmented out code will yields list of densities ->
## subjects -> list of matrixes describing graph
##
## thresholded.matrixces=threshold(netcc, graph.densities)
## thresholded.matrices=llply(graph.densities, function(xx) {
## threshold.correlations(netcc, in.density=xx) },
## .progress=progress.bar.type, .parallel=parallel.executation)

thresholded.matrices=threshold.correlations(netcc, graph.densities)

save.structures.list=append(save.structures.list, "thresholded.matrices")

## ##############################################################################

## load igraph packages
library(igraph)
library(chron)

weighted=FALSE
cat("*** Creating", ifelse(weighted, "WEIGHTED", "UN-WEIGHTED"), "graphs at", length(graph.densities),
    ifelse(length(graph.densities) > 1,
           "different densities", "density"),
    "for each subject\n")
## now create a graph for each subject at each density
if (weighted) {
    g <- llply(thresholded.matrices, lapply,
               function(xx) {
                   graph_from_adjacency_matrix(xx$weights, mode="undirected", weighted=weighted, diag=FALSE)
               },
               .progress=progress.bar.type, .parallel=parallel.executation)
} else {
    g <- llply(thresholded.matrices, lapply,
               function(xx) {
                   graph_from_adjacency_matrix(xx$r.thresh, mode="undirected", diag=FALSE)
               },
               .progress=progress.bar.type, .parallel=parallel.executation)
}

## ##############################################################################
modality="RSFC"
cat("*** Applying set.brainGraph.attributes to each subject's graph at", length(graph.densities), "different densities\n")
## library(brainGraph)

begin.at=Sys.time()
cat(sprintf("***  Starting at %s\n", begin.at))

## g1 <- set.brainGraph.attributes(g[[1]][[1]],
##                                 atlas    = atlas,
##                                 modality = modality,
##                                 subject  = as.character(characteristics.df[1, "Study.ID"]),
##                                 group    = as.character(characteristics.df[1, "Group"]))
## g2 <- set.brainGraph.attributes(g[[length(graph.densities)]][[1]],
##                                 atlas    = atlas,
##                                 modality = modality,
##                                 subject  = as.character(characteristics.df[1, "Study.ID"]),
##                                 group    = as.character(characteristics.df[1, "Group"]))

## this code will only work if g is a list (one element for each
## density) of a list (for element each subject) of graphs (lists)

g.attributes= foreach(dd=icount(length(g))) %do% {
    start=Sys.time()
    cat(sprintf("*** Started subject %s (%03d of %03d) at %s. ",
                names(g)[dd], dd, length(g), start))
    
    ret.graphs=llply(g[[dd]], set.brainGraph.attributes,  ## my.set.attributes,          
          .progress   = progress.bar.type,
          .parallel   = parallel.executation, ## argument to llply
          use.parallel= FALSE, ## to be passed to set.brainGraph.attributes
          atlas       = atlas,
          modality    = modality,
          group       = as.character(characteristics.df[dd, "Group"]),
          subject     = as.character(characteristics.df[dd, "Study.ID"])
          )
    end=Sys.time()
    
    time.taken=format(as.chron(end) - as.chron(start))
    cat (sprintf("Ended at %s. Time taken %s. %0.2f%% completed\n", end, time.taken, (dd/length(g))*100))

    ret.graphs
}
end.at=Sys.time()
total.time.taken=format(as.chron(end.at) - as.chron(begin.at))
cat (sprintf("*** Ended setting brainGraph attributes at %s. Time taken %s\n", end, total.time.taken))

## now make sure that g.attributes has the same names as g for each subject and density
cat("*** Setting names on g.attributes\n")
names(g.attributes) = names(g)
      
## don't set densities as names of list elements since the density of
## the graph may not exactly the same (due to rounding) as the density
## threshold used to threshold the matrices. The buest way to get a
## graph's density is to apply the igraph::graph.density function to a
## graph of get the density attribute (set by
## set.brainGraph.attributes) from the graph itself
##g.attributes = lapply(g.attributes, function(xx) { names(xx) = names(g[[1]]) ; return(xx) })

cat("*** Rounding graph density value to 2 decimal places\n")
g.attributes = lapply(g.attributes, 
    function(ss) {
        lapply(ss,
               function (gg) {
                   gg$density = round(gg$density, 2)
                   return(gg)
               })
    })

actual.densities=lapply(g.attributes, function(xx) {
    vapply(xx, function(yy) { yy$density }, numeric(1), USE.NAMES=FALSE)
})


subject.to.density=do.call(rbind, lapply(seq_along(actual.densities),
    function(ii) {
        data.frame("Study.ID"=
                       rep(names(actual.densities)[ii],
                           length(actual.densities[[ii]])),
                   "density" = actual.densities[[ii]])
    }))

cat("*** Table of density threshold per subject\n")
cat("*** Columns should add to", length(g.attributes), "\n")
cat("*** Rows should add to", length(graph.densities), "\n")
print(addmargins(xtabs(~ Study.ID + density, data=subject.to.density)))

## g.attributes <- Map(
##     function(xx, yy, zz) {
##         llply(xx, set.brainGraph.attributes,  ## my.set.attributes,          
##               .progress=progress.bar.type,
##               .parallel=parallel.executation,
##               use.parallel=FALSE, ## to be passed to set.brainGraph.attributes
##               atlas=atlas,
##               modality=modality, group=yy, subject=zz)
##     },
##     g, as.list(as.character(characteristics.df$Group)), as.list(as.character(characteristics.df$Study.ID)))

save.structures.list=append(save.structures.list, c("g", "g.attributes", "actual.densities", "subject.to.density"))

cat("*** Saving the following data structures to", saved.graph.data.structures.filename, "\n")
cat(paste("+++ ", unlist(save.structures.list), "\n", sep=""), sep="")
pigz.save(list=unlist(save.structures.list), file=saved.graph.data.structures.filename)


## now scale the columns
## graph.feature.columns=grep(roi.label.regexp, colnames(feature.df))
## feature.df[ , graph.feature.columns] = scale( feature.df[ , graph.feature.columns] ) 

## select and reorder the columns
## feature.df=feature.df[,
##     c("subject", "Group", "Gender", "age.in.years", "Verbal", "Performance", "Full", "telomere.final.T.S", "mt.dna",
##       colnames(feature.df)[grep(roi.label.regexp, colnames(feature.df))])]

## feature.df=feature.df[,
##     c("subject", "Group",
##       colnames(feature.df)[grep(roi.label.regexp, colnames(feature.df))])]
## cat("*** There are", dim(feature.df)[1], "subjects\n")

## na.counts=apply(feature.df[, grep("RSFC", colnames(feature.df))], 2, function (xx) { sum(is.na(xx)) } )
## nan.roi.pair.names = convert.feature.vector.labels.to.roi.names(atlas.labels, names(na.counts)[na.counts > 0])
## cat("*** There are", length(nan.roi.pair.names), "ROI pairs with NaN values\n")

## cat("*** There are", length(grep(roi.label.regexp, colnames(feature.df))), "graph-derived features per subject BEFORE removing near zero predictors\n")
## cat("*** Checking for near zero variance predictors\n")

## feature.df = filter.near.zero.variance.predictors(feature.df, in.predictor.tag="RSFC")
## cat("*** There are", length(grep(roi.label.regexp, colnames(feature.df))), "graph-derived features per subject AFTER removing near zero predictors\n")

## feature.df$subject=as.factor(fix.subject.ids(feature.df, "subject"))
## print(convert.feature.vector.labels.to.roi.names(atlas.labels, colnames(feature.matrix)))


