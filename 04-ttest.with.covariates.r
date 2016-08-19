rm(list=ls())
graphics.off()

library(tidyr)

## Reads the seed file and does the (crude) equivalent of BASH variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character(), quiet=TRUE)
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
    name=basename(inSeedPath)
    if (grepl("\\.nii", name)) {
        return(gsub("\\.nii.*", "", name))
    } else if (grepl("\\+tlrc", name)) {
        return(gsub("\\+tlrc.*", "", name))
    } else {
        return (name)
    }
}

readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", ""))
    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
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

list.excluded.subjects <- function(in.df, in.drop.medicated=TRUE) {

    orig.dims=dim(in.df)
    dl=c()
    print.dropped.subject.list <- function (in.drop.list) {
        if (any(in.drop.list)) {
            cat ("--- The following subjects are to dropped from the data frame\n")
            dl=as.vector (  in.df[which(in.drop.list), "ID"])
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

    if (in.drop.medicated) {
        ## ################################################    
        ## Drop medicated subjects (This list was compiled for the TL
        ## paper by Laura. See
        ## rsfcGraphAnalysis/scripts/analyseGroupSubjects.r and
        ## structuralTelomereAnalysis/scripts/analyseStructuralVolumes_Laura.r
        ## Now remove subjects that are NOT medication-naive (Prozac,
        ## Zoloft, Celexa, Citalopram, Klondpin, Seroquil, Cymbalta and
        ## ADHD meds)
        cat("*** Checking for subjects known to have been on medications\n")
        drop.list.medicated=in.df$ID %in%
        c("130", "132", "318", "319", "320",
          "322", "323", "324", "325", "329", "345", "376", "384", "333", "349")
        drop.count=drop.count + sum(drop.list.medicated)
        dropped.subjects.list=c(dropped.subjects.list, print.dropped.subject.list(drop.list.medicated))
        ## in.df = in.df[ ! drop.list, ]
    }
    
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
    drop.list.rum=in.df$ID %in% c("391")
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

add.seed.zscore.files <- function(in.mgd, in.seedName) {

    in.mgd$InputFile=sapply(in.mgd$subject, function (ss) sprintf("%s/%s/rsfc/%s/%s.z-score+tlrc.HEAD", data.dir, ss, seedName, seedName))
    in.mgd$InputFile.exists=vapply(in.mgd$InputFile, file.exists, FALSE)
    return(in.mgd)
}

save.subject.lists <- function (in.zscore.df, in.seed.name) {

    ncl.list.filename=file.path(group.results.dir, sprintf("subject.list.ncl.%s.csv", in.seed.name))
    mdd.list.filename=file.path(group.results.dir, sprintf("subject.list.mdd.%s.csv", in.seed.name))
    combined.list.filename=file.path(group.results.dir, sprintf("subject.list.mddAndNcl.%s.csv", in.seed.name))    

    ncl.list=subset(zscore.df, Group=="NCL", select="subject")
    mdd.list=subset(in.zscore.df, Group=="MDD", select="subject")    
   
    cat("*** Writing NCL subject list to" , ncl.list.filename, "\n")
    ## print(ncl.list)
    write.table(ncl.list, ncl.list.filename, sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    cat("*** Writing MDD subject list to" , mdd.list.filename, "\n")
    ## print(mdd.list)
    write.table(mdd.list, mdd.list.filename, sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)

    cat("*** Writing combined NCL and MDD subject list to" , mdd.list.filename, "\n")
    ## print(zscore.df[, "subject"])
    write.table(zscore.df[, "subject"], combined.list.filename, sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

study.root.dir=file.path(root.dir, "sanDiego/machLearnT1Analysis")
scripts.dir=file.path(study.root.dir, "scripts")
data.dir=file.path(study.root.dir, "data")
group.results.dir=file.path(data.dir, "group.results")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
seeds.data.dir=file.path(data.dir, "seeds")

demographicsFilename=file.path(admin.data.dir, "0-data_entry_current_2014.csv")
demographics=readCsvFile(demographicsFilename)
colnames(demographics)[colnames(demographics)=="Grp"]="Group"
## 378 identified as transgender, was biologically female and not on any hormone therapy
demographics[demographics$ID=="378", "Gender"]="F"
demographics$Gender=droplevels(demographics$Gender)

wasiFilename=file.path(admin.data.dir, "WASI.csv")
wasi=readCsvFile(wasiFilename, inSubjectColumnName="SubID")

seeds.list=lapply(
    sapply(c("juelich_whole_amygdala_seeds.txt",
             "dlpfc_seeds.txt",
             "sgacc_seeds.txt"), function(xx) {file.path(config.data.dir, xx)}),
    readSeedsFile)

## wasi.column.names=c("Verbal", "Performance", "Full")
## only use the WASI Full score as a covariate as perfromance and
## verbal are highly correlated
wasi.column.names=c("Full", "Verbal", "Performance")

covariate.column.names=c("Full")
## covariate.column.names=c("Full", "age.in.years", "Gender")

include.exclude.filename=file.path(data.dir, "include.exclude.0.25.csv")
include.exclude.df=readCsvFile(include.exclude.filename, inSubjectColumnName="subject")
include.exclude.df=include.exclude.df[, c(1, 7,8,9,10,11,12)]
include.exclude.df$subject=as.character(include.exclude.df$subject)
include.exclude.df[include.exclude.df$subject=="300_A", "subject"] = "169/300_A"
include.exclude.df <-  separate(include.exclude.df, subject, into=c("ID", "timepoint"), sep="_", remove=FALSE)

mgd=cbind(include.exclude.df,
    demographics[match(include.exclude.df$ID, demographics$ID), c("Group", "DOB", "MRI", "Gender", "CDRS.tscore")],
    wasi        [match(include.exclude.df$ID, wasi$SubID), wasi.column.names]
          )
colnames(mgd)[length(colnames(mgd))]=wasi.column.names
mgd=fixDates(mgd)
mgd=computeAge(mgd)

mgd=droplevels(subset(mgd, inc.exc=="INCLUDE"))
print(addmargins(with(mgd, table(Group, Gender))))

drop.subject.list=list.excluded.subjects(mgd, TRUE)
mgd=mgd[ ! mgd$ID %in% drop.subject.list, ]
print(addmargins(with(mgd, table(Group, Gender))))

grouping="mddAndNcl"

## new.ttest.arguments=""
## new.ttest.arguments="-resid %s \\\n-CLUSTSIM \\\n-prefix_clustsim resid.CStemp.%s \\\n-tempdir tmp.resid.clustsim.%s"
new.ttest.arguments="\\\n-CLUSTSIM \\\n-prefix_clustsim resid.CStemp.%s"


## seeds=seeds[[1]]
## seed=seeds[1]
## seedName=getSeedName(seed)
## zscore.df=add.seed.zscore.files(mgd, seedName)
## save.subject.lists(zscore.df, seedName)
    
## stop()
for (seeds in seeds.list) {
    for (seed in seeds) {
        seedName=getSeedName(seed)
        
        cat("####################################################################################################\n")
        cat(sprintf("*** Creating covariate file for the %s seed of the %s grouping\n", seedName, grouping))

        zscore.df=add.seed.zscore.files(mgd, seedName)
        
        save.subject.lists(zscore.df, seedName)
        
        ## now center the numeric covariates
        for (col in covariate.column.names) {
            if ( is.numeric(zscore.df[, col]) ) {
                cat("*** Mean centering", col, "\n")            
                zscore.df[, col] = scale(zscore.df[, col], center=TRUE, scale=FALSE)
            } else {
                ##cat(sprintf("*** Skipping centering for %s: Not a numeric column\n", col))
                cat(sprintf("*** WARNING: Converting %s to underlying numeric representation of the R factor and mean centering that\n", col))
                zscore.df[, col] = scale(as.numeric(zscore.df[, col]), center=TRUE, scale=FALSE)            
            }
        }

        rownames(zscore.df)=NULL

        zscore.df=droplevels(zscore.df)

        ## now reorder the columns of zscore.df to InputFile is last and only
        ## pick those columns that we want to correct for in the t-tests
        zscore.df=zscore.df[, c("subject", "Group", covariate.column.names, "InputFile")]

        covariates.filename=file.path(group.results.dir, gsub("..", ".", paste("3dttest.covariates", grouping, seedName, "txt", sep="."), fixed=TRUE))
        cat("*** Writing covariates file to:", covariates.filename, "\n")
        ## the ordering of the columns in the the write command below is
        ## important. It must be subject <covariates> InputFile.

        write.table(zscore.df[, c("subject", covariate.column.names)], file=covariates.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, eol="\n")

        ## cat("*** The data table is as follows:\n")
        ## print(head(zscore.df))

        ## mask=sprintf("../mask.grey.%s.union.masked+tlrc.HEAD", grouping)
        ## mask=sprintf("../mask.grey.%s.union.masked+tlrc.HEAD", "mddAndCtrl")
        mask=file.path(group.results.dir, "/MNI_caez_N27_brain.3mm+tlrc.HEAD")
        prefix=sprintf("ttest.%s.%s.covaried", grouping, seedName)

        setA.group=levels(zscore.df$Group)[1]
        setA.files=paste(zscore.df[zscore.df$Group==setA.group, "InputFile"], collapse=" \\\n")
        setA.labels.and.files = paste(apply(zscore.df[zscore.df$Group==setA.group, c("subject", "InputFile")], 1, paste, collapse=" "), collapse=" \\\n")
        
        setB.group=levels(zscore.df$Group)[2]
        setB.files=paste(zscore.df[zscore.df$Group==setB.group, "InputFile"], collapse=" \\\n")
        setB.labels.and.files = paste(apply(zscore.df[zscore.df$Group==setB.group, c("subject", "InputFile")], 1, paste, collapse=" "), collapse=" \\\n")
        
        three.d.ttest.command = "3dttest++"
        ## use center NONE here as the covariates are already mean
        ## centered before they were written to the covariates file
        if (nchar(new.ttest.arguments) > 0) {
            ## sprintf(new.ttest.arguments, resid.prefix, seedName, seedName),
            ## three.d.ttest.arguments = sprintf("%s \\\n-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",
            three.d.ttest.arguments = sprintf("%s \\\n-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",            
                sprintf("\\\n-CLUSTSIM \\\n-prefix_clustsim CStemp.%s", seedName),            
                mask, prefix,
                setA.group, setA.labels.and.files,
                setB.group, setB.labels.and.files,
                covariates.filename)
        } else {
            three.d.ttest.arguments = sprintf("-mask %s \\\n-prefix %s \\\n-center NONE \\\n-setA %s %s \\\n-setB %s %s \\\n-covariates %s",
                mask, prefix, setA.group, setA.labels.and.files, setB.group, setB.labels.and.files, covariates.filename)
        }
        
        three.d.ttest.command.script.filename=file.path(scripts.dir, gsub("..", ".", sprintf("04-ttest.withCovariates.%s.%s.sh", grouping, seedName), fixed=TRUE))
        cat("*** Writing the 3dttest++ command to:", three.d.ttest.command.script.filename, "\n")

        ## full.three.d.ttest.command=sprintf("cd %s ; mkdir tmp.resid.clustsim.%s ; %s %s", group.results.dir, seedName, three.d.ttest.command, three.d.ttest.arguments)
        full.three.d.ttest.command=sprintf("#!/bin/bash\n\nunset AFNI_COMPRESSOR \n\nexport OMP_NUM_THREADS=40\n\nmkdir -p %s/ttest.%s \n\ncd %s/ttest.%s \n\n%s %s",
            group.results.dir, seedName,
            group.results.dir, seedName,        
            three.d.ttest.command, three.d.ttest.arguments)    
        cat (full.three.d.ttest.command, file=three.d.ttest.command.script.filename)
        ## cat(full.three.d.ttest.command, "\n\n")
        
    }
}
