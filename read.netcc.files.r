rm(list=ls())
graphics.off()

fix.subject.ids <- function (in.subject.ids) {

    in.subject.ids=gsub("_[ABCD]", "", as.character(in.subject.ids), fixed=FALSE)
    subjectIds=as.character(in.subject.ids)
    if (length(pmatch(subjectIds, "300")) > 0 ) {
        subjectIds=sub("^300", "169/300", subjectIds, fixed=FALSE)
        subjectIds=sub("^0", "", subjectIds, fixed=FALSE)
    }
    return (subjectIds)
}


read.demographics.table <- function(in.filename) {
    cat("*** Reading", in.filename, "\n")
    demographics=read.csv(in.filename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
    cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
    demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))
    demographics=demographics[grep("Test_Subject", demographics$ID, invert=TRUE), ]
    demographics$ID=droplevels(demographics$ID)

    return(demographics)
}


read.wasi.table <- function(in.filename) {
    cat("*** Reading", in.filename, "\n")
    wasi=read.csv(in.filename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
    cat(sprintf("*** Read WASI-II data for %s unique subjects\n",  length(unique(wasi$SubID))))

    return(wasi)
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


filter.nonexistant.files <- function (in.file.list, in.print.nonexistant=FALSE) {
    existant=sapply(in.file.list, file.exists)
    if (in.print.nonexistant) {
        cat("*** The following netcc files do not exist:\n")
        cat(paste("---", in.file.list[!existant]), sep="\n")
    }
    
    aa=in.file.list[existant]

    return(list("filenames"=aa, "subjectNames"=names(aa)))
}

read.netcc.files <- function (in.netcc.filenames, in.subject.names, in.label.names) {
    dnames=list()
    dnames[[1]] = in.label.names$FvLabel
    dnames[[2]] = in.label.names$FvLabel
    dnames[[3]] = in.subject.names

    netcc=array(as.numeric(NA),
        dim=c(dim(in.label.names)[1], dim(in.label.names)[1], length(in.subject.names)),
        dimnames=dnames)
    ## cat("*** netcc array is", paste(dim(netcc)), "\n")
    for (ii in seq.int(1, length(in.subject.names))) {
        cat(sprintf("+++ Reading (%04d/%04d) %s\r", ii, length(in.subject.names), in.netcc.filenames[ii]))
        ff=read.table(in.netcc.filenames[ii], header=FALSE, skip=6)
        netcc[, , ii ] = as.matrix(ff)
    }
    cat("\n")
    return(netcc)
}

read.aal.atlas.label.names <- function (in.filename) {
    cat("*** Reading", in.filename, "\n")
    aal.labels=read.table(in.filename, header=FALSE)
    ##    aal.labels[, 3] = apply(aal.labels, 1, function(xx) { sprintf("R%03d.%s", as.integer(xx[1]), xx[2]) })
    aal.labels[, 3] = apply(aal.labels, 1, function(xx) { sprintf("R%03d", as.integer(xx[1])) })    
    colnames(aal.labels) = c("ID", "Name", "FvLabel")

    return(aal.labels)
}

convert.to.feature.matrix <- function(in.netcc, in.correlation.tag) {

    ut=upper.tri(in.netcc[, , 1])
    
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
    ## the labels will be of the form R000.R000 whcih can be
    ## translated back to ROI names by using the integer component to
    ## index the aal labels data frame
    colnames(fv)=apply(nn, 1, function(xx) { sprintf("%s.%s.%s", in.correlation.tag, rn[xx[1]], cn[xx[2]]) } )
    
    return(fv)
}

convert.feature.vector.labels.to.roi.names <- function (in.aal.labels, in.feature.vector.names) {

    sapply(in.feature.vector.names,
           function(xx) {
               indices=as.integer(strsplit(gsub(roi.label.regexp, "\\2 \\3", xx), " ")[[1]])
               paste(in.aal.labels[indices, 2], collapse=" <-> ")
           })

}

filter.subjects <- function(in.df) {

    print.dropped.subject.list <- function (in.drop.list) {
        if (any(in.drop.list)) {
            cat ("--- The following subjects were dropped from the data frame\n")
            cat ("---", paste (as.vector (  in.df[which(in.drop.list), "subject"]), collapse=" "), "\n")
            cat ("---", paste (as.vector (  in.df[which(in.drop.list), "Grp"]),     collapse=" "), "\n")
        }    
    }
    
    cat("*** Subject numbers BEFORE filtering\n")
    print(addmargins(table(in.df[, c("Gender", "Grp")])))

    ## ################################################
    ## Drop MDDs with CDRS-R < 55
    cat("*** Checking for NCLs with CDRS-R tscore > 54\n")
    drop.list=in.df$Grp=="NCL" & in.df$CDRS.tscore > 54
    print.dropped.subject.list(drop.list)
    in.df = in.df[ ! drop.list , ]

    ## ################################################    
    ## Drop NCLs with CDRS-R > 55
    cat("*** Checking for MDDs with CDRS-R tscore < 55\n")
    drop.list=in.df$Grp=="MDD" & in.df$CDRS.tscore < 55
    print.dropped.subject.list(drop.list)
    in.df = in.df[ ! drop.list, ]    

    ## ################################################    
    ## Drop medicated subjects (This list was compiled for the TL
    ## paper. See rsfcGraphAnalysis/scripts/analyseGroupSubjects.r and
    ## structuralTelomereAnalysis/scripts/analyseStructuralVolumes_Laura.r
    ## Now remove subjects that are NOT medication-naive (Prozac, Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbalta and ADHD meds)
    drop.list=in.df$subject %in% c("130", "132", "318", "319", "320", "322", "323", "324", "325", "329", "345", "376", "384", "333", "349")
    print.dropped.subject.list(drop.list)
    in.df = in.df[ ! drop.list, ]

    ## ################################################    
    ## Drop subjects with WASI Performance and WASI Verbal < 80
    cat("*** Checking for WASI Verbal < 80 and WASI Perfromance < 80\n")
    drop.list=in.df$Verbal < 80 & in.df$Performance < 80
    print.dropped.subject.list(drop.list)
    in.df = in.df[ ! drop.list, ]

    cat("*** Checking for other excluded subjects\n")
    ## 391 is in group RUM on the demographics sheet, but I've no idea
    ## what this means so dumo them
    drop.list=in.df$subject %in% c("391")
    print.dropped.subject.list(drop.list)
    in.df = in.df[ ! drop.list, ]
    
    cat("*** Subject numbers AFTER filtering\n")
    print(addmargins(table(in.df[, c("Gender", "Grp")])))
    
    return(in.df)
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

study.root.dir=file.path(root.dir, "sanDiego/machLearnT1Analysis")
standard.data.dir=file.path(study.root.dir, "standard")

scripts.dir=file.path(study.root.dir, "scripts")
data.dir=file.path(study.root.dir, "data")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")

## regular expression to match the labels given to the elements in the
## feature vector of each subject
roi.label.regexp="(RSFC|GM|DTI)\\.R([0-9]{3})\\.R([0-9]{3})"

subjects=dir(data.dir, pattern="[0-9][0-9][0-9]_A*")

demographics.filename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=read.demographics.table(demographics.filename)

## 378 identified as transgender, was biologically female and not on any hormone therapy
demographics[demographics$ID=="378", "Gender"]="F"
demographics$Gender=droplevels(demographics$Gender)

atlas.filename=file.path(standard.data.dir, "aal2_for_SPM12", "aal2.nii.txt")
atlas.labels=read.aal.atlas.label.names(atlas.filename)

telomere.filename=file.path(admin.data.dir, "Tony Yang TL data 06162014.csv")
telomere.data=read.telomere.table(telomere.filename)

mt.dna.filename=file.path(admin.data.dir, "mtDNA.csv")
mt.dna.data=read.mt.dna.table(mt.dna.filename)

wasi.filename=file.path(admin.data.dir, "WASI.csv")
wasi.data=read.wasi.table(wasi.filename)


## now read the graphs created from teh RSFC analysis (note the
## "rsfcGraphs" argument)
netcc.filenames.and.subejcts=filter.nonexistant.files(build.netcc.filenames(subjects, "rsfcGraphs"), in.print.nonexistant=TRUE)
## print(netcc.filenames.and.subejcts$filenames)
## print(netcc.filenames.and.subejcts$subjectNames)

netcc=read.netcc.files(netcc.filenames.and.subejcts$filenames, netcc.filenames.and.subejcts$subjectNames, atlas.labels)

## feature.matrix=convert.to.feature.matrix(netcc[1:5, 1:5, 1:5], "RSFC")
feature.matrix=convert.to.feature.matrix(netcc, "RSFC")
## subject.list=fix.subject.ids(netcc.filenames.and.subejcts$subjectNames[1:5])
subject.list=fix.subject.ids(netcc.filenames.and.subejcts$subjectNames)
feature.df=data.frame(
    "subject"           =as.factor(subject.list),
    demographics                      [match(subject.list, demographics$ID),         c("Grp", "Gender", "DOB", "MRI", "CDRS.tscore")],
    wasi.data                         [match(subject.list, wasi.data$SubID),         c("Verbal", "Performance", "Full")],
    "telomere.final.T.S"=telomere.data[match(subject.list, telomere.data$sample.ID), "final.T.S"],
    "mt.dna"            =mt.dna.data  [match(subject.list, mt.dna.data$subject),     "mtDNA"], 
    feature.matrix)

feature.df=fix.dates(feature.df)
feature.df=compute.age(feature.df)

## remove medicated subjects, subjects with too high/low CDRS-R and WASI score
feature.df=droplevels(filter.subjects(feature.df))

## reorder the columns
feature.df=feature.df[,
    c("subject", "Grp", "Gender", "DOB", "MRI", "age.in.years", "Verbal", "Performance", "Full", "telomere.final.T.S", "mt.dna",
    colnames(feature.df)[grep(roi.label.regexp, colnames(feature.df))])]

## feature.df$subject=as.factor(fix.subject.ids(feature.df, "subject"))
## print(convert.feature.vector.labels.to.roi.names(atlas.labels, colnames(feature.matrix)))


