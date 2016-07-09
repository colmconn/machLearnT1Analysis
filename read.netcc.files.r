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
        cat("+++", str_wrap(paste (names( zero.subjects)[ zero.subjects == FALSE ], collapse=" "), width=80), "\n")

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
            cat ("--- The following subjects were dropped from the data frame\n")
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

## ##################################################################################################
## END OF FUNCITONS
## ##################################################################################################

### setup path variables
if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

library(stringr)
library(tidyr)
library(igraph)
library(plyr)
library(data.table)
library(brainGraph)

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

atlas="brainGraph::aal2_94"
atlas.dt=eval(parse(text=atlas))
cat("*** Using the", atlas, "atlas\n")

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
    save(netcc, netcc.filenames.and.subjects.df, file=saved.netcc.filename)
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

## drop.subject.list=paste(drop.subject.list, "_A", sep="")
## remove subjects from the netcc array
## cat("*** BEFORE removing subjects from netcc", dim(netcc), '\n')
netcc=drop.subjects.from.netcc(netcc, drop.subject.list)
## print(dim(netcc))
## cat("*** AFTER removing subjects from netcc", dim(netcc), '\n')

## remove subjects from the netcc.filenames.and.subjects.df
## cat("*** BEFORE removing subjects from netcc.filenames.and.subjects.df", dim(netcc.filenames.and.subjects.df), '\n')
netcc.filenames.and.subjects.df=drop.subjects(netcc.filenames.and.subjects.df, drop.subject.list)
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

## now that all of the filtering has been done it's time to z-score the array
cat("*** Z-scoring netcc array\n")
netcc.og=netcc
netcc=atanh(netcc)

cat("*** subject distribution is as follows:\n")
subject.distribution=addmargins(xtabs(~ Group + Gender, data=characteristics.df))
print(subject.distribution)

threshold.correlations <- function (in.netcc, in.density=0.1) {
    subject.count=dim(in.netcc)[3]
    N <- ncol(in.netcc[, , 1])
    emax <- N  * (N - 1) / 2

    out=list()
    for (ii in seq.int(1, subject.count)) {

        ## the following code is from corr_matrix.R in brainGraph
        ##
        ## try to make the data structures as similar as possible, if
        ## not identical, so that we can use brainGraph for other
        ## processess later if necessary
        ## r=in.netcc[ , , ii]
        thresh=sort(in.netcc[ , , ii][lower.tri(in.netcc[ , , ii])])[emax - in.density * emax]
        r.thresh <- ifelse(in.netcc[ , , ii] > thresh, 1, 0)
        out[[ii]] <- list(R=in.netcc[ , , ii], P=pnorm(netcc[ , , ii]), r.thresh=r.thresh, threshold=thresh)
    }
    names(out)=names(in.netcc[1, 1, ])
    return(out)
}
cat("*** Thresholding correlations for each subject\n")
thresholded.matrixes=threshold.correlations(netcc)

cat("*** Creating graph for each subject\n")
## now create a graph for each subject
g <- lapply(thresholded.matrixes, function(xx) { graph_from_adjacency_matrix(xx$r.thresh, mode="undirected", diag=FALSE) } )


## modality="RSFC"

## g1 <- set.brainGraph.attributes(g[[1]],
##                                 modality=modality,
##                                 subject=characteristics.df[1, "Study.ID"],
##                                 group=characteristics.df[1, "Group"])

## g <- Map(function(x, y, z) {
##     llply(x, set.brainGraph.attributes, atlas=atlas, modality=modality, group=y, subject=z, .progress='text')
## },
##          g, as.list(as.character(characteristics.df$Group)), as.list(as.character(characteristics.df$Study.ID)))



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


