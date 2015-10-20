rm(list=ls())

library(reshape)

graphics.off()

isExcluded <- function(sid) {
    ## cat(sid, "\n")
    if ( sid %in% c("311_A", "119_A","346_A", "370_A", "130_A", "132_A", "319_A", "320_A", "322_A", "323_A", "325_A", "329_A", "333_A", "376_A", "378_A") )
        return (TRUE)
    else
        return (FALSE)
}

getGroup <- function(sid) {
    if (sid == "300")
        sid="169/300"
    
    return(as.character(demographics[match(sid, demographics$ID), "Grp"]))
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

makeDesignMatrix <- function(inData, inEvs, inCovariates=NULL) {

    checkVarsInDataFrame <- function(inData, inNames) {
        return(all(inNames %in% colnames(inData)))
    }

    if ( is.null(inEvs) ) {
        stop("You did not provide the names of any explanatory variables. Cannot continue.\n")
    }
        
    if ( ! is.null(inEvs) & ! checkVarsInDataFrame( inData, inEvs ) ) {
        stop("Explanatory variable names are not in the column names of the data frame. Cannont continue.\n")
    }
    
    if ( ! is.null(inCovariates) & ! checkVarsInDataFrame( inData, inEvs ) ) {
        stop("Explanatory variable names are not in the column names of the data frame. Cannont continue.\n")
    }

    mx=as.matrix(inData[c(inEvs, inCovariates)])
    rownames(mx)=NULL

    return(mx)
}


demean <- function (inData, inColumnNames, verbose=FALSE) {

    for ( cname in inColumnNames ) {
        if ( cname %in% colnames(inData) ) {
            mn=mean(inData[, cname])
            if (verbose) {
                cat(paste("*** Mean for", cname, "is", mn, "\n"))
            }
            inData[, cname] = inData[, cname] - mn
        } else {
            cat(paste(cname, "not in the columns of data frame. Skipping\n"))
        }
    }
    return(inData)
}


write.model.matrix <- function (inData) {
    write.table(inData, stdout(), quote=FALSE, col.names=TRUE, row.names=FALSE)
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

admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/admin"))
config.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/config"))
group.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.data"))
group.results.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/Group.results"))
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts"))
##vbm.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm/struc"))

## vbm.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbmFromFreesurfer/struc"))

vbm.data.dir=normalizePath(file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/vbm.subject.list.from.matthew.n111/struc"))

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ""))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

wasiFilename=file.path(admin.data.dir, "WASI.csv")
cat("*** Reading", wasiFilename, "\n")
wasi=read.csv(wasiFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read WASI-II data for %s unique subjects\n",  length(unique(wasi$SubID))))


tiv.filename=file.path(vbm.data.dir, "tiv.csv")
if ( file.exists ( tiv.filename ) ) {
    tiv=read.csv(tiv.filename, header=TRUE)
    tiv$subject=gsub("((?:MDD|NCL)[.])?(?:([0-9]+)_A[0-9]?)", "\\2", as.character(tiv$subject), fixed=FALSE)
    tiv$subject=as.factor(gsub("300", "169/300", as.character(tiv$subject), fixed=TRUE))
} else {
    cat("### Run calculateTiv.sh first to calculate TIV and TBV if you want them in the design matrixes.\n")
}

createSubjectListFiles=FALSE
renameFiles=FALSE

## subjects=read.table(file.path(config.data.dir, "subject_list_from_matthew.txt"), header=FALSE)

subjects = dir(vbm.data.dir, pattern=".*[_.]anat_struc.nii.gz")
if (length(subjects) == 0 ) {
    stop(paste("No files in", vbm.data.dir, "matched the file selection pattern. Stopping.\n"))
}

subjects=data.frame(
    "filename"=subjects,
    "subject"=gsub("((?:MDD|NCL)[.])?(?:([0-9]+)_A[0-9]?)[_.]anat_struc.nii.gz", "\\2", subjects, fixed=FALSE))
rownames(subjects)=NULL
subjects$subject=as.factor(gsub("300", "169/300", as.character(subjects$subject), fixed=TRUE))


if (exists("tiv")) {
    subjects=cbind(
        subjects,
        demographics[match(subjects$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "CDRS.tscore")],
        wasi        [match(subjects$subject, wasi$SubID), c("Verbal", "Performance", "Full")],
        tiv         [match(subjects$subject, tiv$subject), c("tiv", "tbv")])
    
} else {
    subjects=cbind(
        subjects,
        demographics[match(subjects$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "CDRS.tscore")],
        wasi        [match(subjects$subject, wasi$SubID), c("Verbal", "Performance", "Full")])
}
subjects <- rename(subjects, c("Grp"="group"))

subjects$NCL=ifelse(subjects$group=="NCL", 1, 0)
subjects$MDD=ifelse(subjects$group=="MDD", 1, 0)

subjects$Gender.num=ifelse(subjects$Gender=="F", 1, 0)

subjects=fixDates(subjects)
subjects=computeAge(subjects)
rownames(subjects)=NULL

## now demean the gender and age.in.years columns
if (exists("tiv")) {
    subjects.demeaned=demean(subjects, c("age.in.years", "Gender.num", "CDRS.tscore", "Full", "tiv", "tbv"), verbose=TRUE)
} else {
    subjects.demeaned=demean(subjects, c("age.in.years", "Gender.num", "CDRS.tscore", "Full"), verbose=TRUE)
}

##################################################
### setup the variables with the design matrixes
##################################################
simple.two.group.mx=makeDesignMatrix(subjects, c("MDD", "NCL"))
simple.two.group.mx.with.covariates=makeDesignMatrix(subjects.demeaned, c("MDD", "NCL"), c("age.in.years", "Gender.num", "Full"))
if (exists("tiv")) {
    simple.two.group.mx.with.wasi.and.tbv=makeDesignMatrix(subjects.demeaned, c("MDD", "NCL"), c("Full", "tbv"))
} else {
    simple.two.group.mx.with.wasi=makeDesignMatrix(subjects.demeaned, c("MDD", "NCL"), c("Full"))
}

## now set up the MDD-only regression matrices
mdd.subjects.only=subset(subjects, group=="MDD")
mdd.subjects.only.demeaned=demean(mdd.subjects.only, c("age.in.years", "Gender.num"), verbose=TRUE)
cdrsr.regression.mx.withcovariates=makeDesignMatrix(mdd.subjects.only.demeaned, c("CDRS.tscore"), c("age.in.years", "Gender.num"))

##################################################

if (createSubjectListFiles) {
    mdd.subjects.only$subject=sapply(mdd.subjects.only$subject,
        function (xx) {
            ifelse(xx=="117" | xx=="334",
                   paste(xx, "A2", sep="_"), paste(xx, "A", sep="_"))
        })
    mdd.subjects.only$subject=gsub("169/300", "300", as.character(mdd.subjects.only$subject), fixed=TRUE)
    write.table(mdd.subjects.only$subject, file.path(vbm.data.dir, "../mdd.subjectlist.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)

    ncl.subjects=subset(subjects, group=="NCL", select=subject)
    ncl.subjects$subject=sapply(ncl.subjects$subject,
        function (xx) {
            ifelse(xx=="117" | xx=="334",
                   paste(xx, "A2", sep="_"), paste(xx, "A", sep="_"))
        })
    ncl.subjects$subject=gsub("169/300", "300", as.character(ncl.subjects$subject), fixed=TRUE)
    write.table(ncl.subjects$subject, file.path(vbm.data.dir, "../ncl.subjectlist.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
}


if (renameFiles) {                     

    owd=getwd()
    setwd(vbm.data.dir)
    files.to.be.renamed=c(dir(".", pattern="^[0-9]+_A[0-9]?[_.]anat_struc.nii.gz"), dir(".", pattern="^[0-9]+_A[0-9]?[_.]anat_struc_brain.nii.gz"))
    if (length(files.to.be.renamed) > 0 )  {
        for (ff in files.to.be.renamed ) {
            subject=gsub("_A[0-9]?[_.]anat_struc.*", "", ff, fixed=FALSE)
            group=getGroup(subject)
            ## cat(subject, "\n")
            
            if(isExcluded(subject)) {
                warning(subject, "is excluded. Deleting", ff)
                ## file.remove(ff)
            } else if (group == "MDD" | group == "NCL") {
                newname=paste(group, ff, sep=".")
                cat("Renaming", ff, "to", newname, "\n")
                file.rename(ff, newname)
            } else {
                cat("*** Don't know what to do with", subject, "***\n")
            }
        }
    } else {
        cat(paste("*** No files in", getwd(), "to be renamed\n"))
    }
    setwd(owd)

    owd=getwd()
    setwd(file.path(vbm.data.dir, "../"))
    files.to.be.renamed=c(dir(".", pattern="^[0-9]+_A[0-9]?[_.]anat.nii.gz"))
    if (length(files.to.be.renamed) > 0) {
        for (ff in files.to.be.renamed ) {
            subject=gsub("_A[0-9]?[_.]anat.*", "", ff, fixed=FALSE)
            group=getGroup(subject)
            ## cat(subject, "\n")
            
            if(isExcluded(subject)) {
                warning(subject, "is excluded. Deleting", ff)
                ## file.remove(ff)
            } else if (group == "MDD" | group == "NCL") {
                newname=paste(group, ff, sep=".")
                cat("Renaming", ff, "to", newname, "\n")
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
