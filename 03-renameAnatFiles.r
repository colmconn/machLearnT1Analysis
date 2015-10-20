#!/bin/env Rscript

library(reshape)
library(plyr)

isMedicated <- function(inData) {
    ## cat(sid, "\n")

    ## the following medicaton-based exclusion criteria are taken from Laura's /data/sanDiego/structuralTelomereAnalysis/scripts/n=107.trending.TL.Volume.r
    ## Now remove subjects that are NOT medication-naive (Prozac, Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbatta) (n=132)/(n=125)
    #mgd <- mgd[! mgd$subject %in% c("130", "132", "318", "319", "320", "322", "323", "324", "325", "329", "345", "376", "384"), ]
    ## Now remove subjects that were on ADHD medication (Focolin, Stratera) (n=125)/(n=123)
    #mgd <- mgd[! mgd$subject %in% c("333", "349"), ]

    return ( inData$subject %in% c("119", "130", "132", "311", "319", "320", "322", "323", "325", "329", "333", "346", "370", "376", "378"))
}

isControlCdrsrTscoreGreaterThanFiftyFour <- function (inData) {
    return (inData$group == "NCL" & inData$CDRS.tscore > 54)
}

isWasiLessThanEighty <- function(inData) {
    return (inData$Verbal < 80 & inData$Performance < 80)
}    
    
determineIfIncluded <- function (inData) {


    inData$medicated=isMedicated(inData)
    inData$ncl.cdrsr.score.gt.54 = isControlCdrsrTscoreGreaterThanFiftyFour(inData)
    inData$wasiTooLow=isWasiLessThanEighty(inData)

    ## ordinarily you'd use apply here instead of alply, however using
    ## apply coerses the data.frame to a matrix and in so doing calls
    ## format which for some inconceivable reason preprends a space to
    ## the string representation of TRUE whcih is not recognized as
    ## TRUE when coersion back to a logical takes place. Try
    ## format(c(TRUE, FALSE)) at the R prompt to see. Also see
    ## http://stackoverflow.com/questions/18614236/apply-prepends-space-for-logical
    inData$include = unlist(alply(inData, 1, function(xx) {  ! any(unlist(xx[c("medicated", "ncl.cdrsr.score.gt.54", "wasiTooLow")]))  } ))

    return (inData)
}

fixDates <- function (inData) {
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

computeAge <- function(inData) {

    age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inData$age.in.years=age.in.weeks/52

    return(inData)
}

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

makeListOfAnatFiles <- function() {
    ## 157_A.anat_struc.nii.gz
    subjects = dir(vbm.data.dir, pattern=".*[_.]anat_struc.nii.gz")
    if (length(subjects) == 0 ) {
        stop(paste("No files in", vbm.data.dir, "matched the file selection pattern. Stopping.\n"))
    }
    return (subjects)
}

getGroup <- function(sid) {
    if (sid == "300")
        sid="169/300"
    
    return(as.character(subjects[match(sid, subjects$subject), "group"]))
}


####################################################################################################
### Main code
####################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

study.root.dir=normalizePath(file.path(root.dir, "sanDiego/machLearnT1Analysis"))
admin.data.dir=normalizePath(file.path(study.root.dir, "data/admin"))
## config.data.dir=normalizePath(file.path(study.root.dir, "data/config"))
scripts.dir=normalizePath(file.path(study.root.dir, "scripts"))

vbm.root.dir=normalizePath(file.path(study.root.dir, "data/vbm"))
vbm.data.dir=normalizePath(file.path(vbm.root.dir, "struc"))

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)

subjects=makeListOfAnatFiles()
wasiFilename=file.path(admin.data.dir, "WASI.csv")
wasi=readCsvFile(wasiFilename, "SubID")

subjects=data.frame(
    "filename"=subjects,
    "subject"=gsub("((?:MDD|NCL)[.])?(?:([0-9]+)_A[0-9]?)[_.]anat_struc.nii.gz", "\\2", subjects, fixed=FALSE),
    "ID" = gsub("((?:MDD|NCL)[.])?(([0-9]+)_A[0-9]?)[_.]anat_struc.nii.gz", "\\2", subjects, fixed=FALSE),
    "timepoint"= gsub("((?:MDD|NCL)[.])?(([0-9]+)_(A)[0-9]?)[_.]anat_struc.nii.gz", "\\4", subjects, fixed=FALSE))
rownames(subjects)=NULL
subjects$subject=as.factor(gsub("300", "169/300", as.character(subjects$subject), fixed=TRUE))

subjects=cbind(
    subjects,
    demographics[match(subjects$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "CDRS.tscore")],
    wasi        [match(subjects$subject, wasi$SubID), c("Verbal", "Performance", "Full")])
subjects=rename(subjects, c("Grp"="group"))
subjects=computeAge(fixDates(subjects))
rownames(subjects) = NULL
subjects=determineIfIncluded(subjects)


if (any(subjects$group=="RUM") ) {
    ## We'll assume that the subject (391) with the group name RUM is
    ## an MDD since their CDRS-R score is 76
    subjects[subjects$group=="RUM", "group"]="MDD"
}

##print (subjects[is.na(subjects$include) | subjects$include==FALSE, ])
##print(dim(subjects[is.na(subjects$include) | subjects$include==FALSE, ]))
##print (subjects)

####################################################################################################
### move files belonging to excluded subjects
####################################################################################################
excluded=subset(subjects, include==FALSE | is.na(include))
excluded.path=file.path(vbm.root.dir, "excluded")
struc.excluded.path=file.path(excluded.path, "struc")
if ( ! dir.exists(excluded.path) ) {
    dir.create(excluded.path)
    dir.create(struc.excluded.path)
}

files.to.be.moved[["from"]]=c(
                     sapply(excluded$ID, function(xx) { file.path(vbm.root.dir, paste(xx, "anat.nii.gz",             sep=".")) }),
                     sapply(excluded$ID, function(xx) { file.path(vbm.data.dir, paste(xx, "anat_struc.nii.gz",       sep=".") ) } ),
                     sapply(excluded$ID, function(xx) { file.path(vbm.data.dir, paste(xx, "anat_struc_brain.nii.gz", sep=".") ) } ) )

files.to.be.moved[["to"]]  =c(
                     sapply(excluded$ID, function(xx) { file.path(excluded.path,       paste(xx, "anat.nii.gz",             sep=".")) }),
                     sapply(excluded$ID, function(xx) { file.path(struc.excluded.path, paste(xx, "anat_struc.nii.gz",       sep=".") ) } ), 
                     sapply(excluded$ID, function(xx) { file.path(struc.excluded.path, paste(xx, "anat_struc_brain.nii.gz", sep=".") ) } ) )

if (length(files.to.be.moved[["from"]]) > 0 ) {
    cat("*** The following files belonging to excluded subjects are being moved to the excluded directory:\n")
    for (ii in 1:length(files.to.be.moved[["from"]]) ) {
        cat(sprintf("%s => %s\n", files.to.be.moved[["from"]][ii], files.to.be.moved[["to"]][ii]))
    }
    file.rename(files.to.be.moved[["from"]], files.to.be.moved[["to"]])
} else {
    cat("*** No files belonging to excluded subjects found to be moved\n")
}

####################################################################################################
## rename the files
####################################################################################################
directory.to.suffix.mapping=list(
    list("directory"=vbm.root.dir, "suffix"="anat.nii.gz"),
    list("directory"=vbm.data.dir, "suffix"="anat_struc.nii.gz"),
    list("directory"=vbm.data.dir, "suffix"="anat_struc_brain.nii.gz"))

for (ii in 1:length(directory.to.suffix.mapping)) {
    directory=directory.to.suffix.mapping[[ii]]$directory
    suffix=directory.to.suffix.mapping[[ii]]$suffix

    owd=getwd()
    setwd(directory)

    files.to.be.renamed=dir(".", pattern=paste("^[0-9]+_A[0-9]?[_.]", suffix, sep=""))

    if (length(files.to.be.renamed) > 0 )  {
        cat("*** Renaming files in", directory, "\n")
        for (ff in files.to.be.renamed ) {
            subject=gsub(paste("_A[0-9]?[_.]", suffix, sep=""), "", ff, fixed=FALSE)
            group=getGroup(subject)
            
            if (group == "MDD" || group == "NCL") {
                newname=paste(group, ff, sep=".")
                cat("*** Renaming", ff, "=>", newname, "\n")
                file.rename(ff, newname)
            } else {
                cat("*** Don't know what to do with", subject, "***\n")
            }
        }
    } else {
        cat(paste("*** No files in", getwd(), "to be renamed\n"))
    }
    setwd(owd)
}        


####################################################################################################
### make the subject list files
####################################################################################################
mdd.subjects=gsub("[A-Z]{3}[._]([0-9]+_A[0-9]?).*", "\\1", dir(vbm.root.dir, "^MDD.*"), fixed=FALSE)
ncl.subjects=gsub("[A-Z]{3}[._]([0-9]+_A[0-9]?).*", "\\1", dir(vbm.root.dir, "^NCL.*"), fixed=FALSE)

if (length(mdd.subjects) > 0) {
    subject.list.filename=file.path(vbm.root.dir, "mdd.subjectlist.txt")
    cat("*** Writing", subject.list.filename, "\n")
    write.table(mdd.subjects, subject.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
}
if (length(ncl.subjects) > 0) {
    subject.list.filename=file.path(vbm.root.dir, "ncl.subjectlist.txt")
    cat("*** Writing", subject.list.filename, "\n")        
    write.table(ncl.subjects, subject.list.filename, quote=FALSE, col.names=FALSE, row.names=FALSE)
}


####################################################################################################
### create the template list file
####################################################################################################
template.list.command=paste("./createTemplateList.pl", vbm.root.dir)
cat("*** Creating the template list file using the command:", template.list.command, "\n")
system(template.list.command)

####################################################################################################
### make simple two group design matrix and contrast for use with the
### very first invokation of randomise in fslvbm_3_proc
####################################################################################################
n.mdd.subjects=length(mdd.subjects)
n.ncl.subjects=length(ncl.subjects)
design.command=sprintf("design_ttest2 design %d %d", n.mdd.subjects, n.ncl.subjects)
cat("*** Creating simple two group design matrixc and contrast using the command:", design.command, "\n")
owd=getwd()
setwd(vbm.root.dir)
system(design.command)
setwd(owd)

####################################################################################################
### Inform on how to run the rest of the VBM
####################################################################################################
cat("*** You can now run the fslvbm_2_template -n in", vbm.root.dir, "\n")
cat("*** When that command finishes and you inspect the template, you can then run the fslvbm_3_proc in", vbm.root.dir, "\n")
