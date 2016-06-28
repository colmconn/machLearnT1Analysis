#!/bin/env Rscript

options(error = function() traceback())

AFNI_R_DIR=Sys.getenv("AFNI_R_DIR", unset=NA)

## use the functions for loading and saving briks from the AFNI
## distribution as they can cleanly handle floats/doubles
if ( ! is.na(AFNI_R_DIR) ) {
    source(file.path(AFNI_R_DIR, "AFNIio.R"))
} else {
    stop("Couldn't find AFNI_R_DIR in environment. This points to the location from which to load functions for reading and writing AFNI BRIKS. Stopping!")
}

library(getopt)
##########################################################################################################################################################################
### START OF FUNCTIONS ###################################################################################################################################################
##########################################################################################################################################################################


help <- function(){

}

check.command.line.arguments <- function (in.opt) {

    if (is.null(in.opt$quiet)) {
        in.opt$quiet=FALSE
    }

    if (is.null(in.opt$execute)) {
        in.opt$execute=FALSE
    }
    
    ## if help was asked for print a friendly message
    ## and exit with a non-zero error code
    if ( !is.null(in.opt$help) ) {
        cat(getopt(spec, usage=TRUE));
        if ( ! interactive()) 
            q(status=1)
    }

    if ( is.null(in.opt$window) ) {
        if (! in.opt$quiet)
            cat("*** No value provided for window. Assuming none\n")
        in.opt$window="none"
    }
    
    valid.window.types=c("consecutive", "overlap", "none")
    window.type=pmatch(in.opt$window, valid.window.types)
    if ( is.na(window.type) ) {
        cat(sprintf("*** Valid options for window type are %s\n", paste(valid.window.types, collapse=", ")))
        if ( ! interactive()) 
            q(status=1)
    } else {
        in.opt$window=valid.window.types[window.type]
    }

    if (in.opt$window == "none" ) 
        in.opt$width=NA
    
    if ( is.null(in.opt$width) && in.opt$window != "none") {
        cat("*** You must provide a window width when using a \"", in.opt$window,  "\" window\n", sep="")
        if ( ! interactive()) 
            q(status=1)
        
    }##  else {
        ## in.opt$width=NA
    ## }

    if ( is.null(in.opt$step) && in.opt$window == "overlap") {
        cat("*** You must provide a step to leave between consecutive window starts when using a \"", in.opt$window,  "\" window\n", sep="")
        q(status=1)
    } else if ( in.opt$window %in% c("consecutive", "none") ) {
        if (! in.opt$quiet)
            cat("*** Window type is set to \"", in.opt$window, "\", ignoring step\n", sep="")
        in.opt$step = NA
    }

    if ( is.null(in.opt$source)) {
        cat("*** You must provide a source file\n")
        if ( ! interactive()) 
            q(status=1)
    }

    if (is.null(in.opt$destination)) {
        cat("*** You must provide a destination directory. It will be created if it does not already exist.\n")
        if ( ! interactive()) 
            q(status=1)
    }

    if (! is.null(in.opt$extra)) {
        in.opt$extra=gsub("[\\\'\"]", "", in.opt$extra)
    }


    
    return(in.opt)
}

print.command.line.arguments.summary <- function () {
    cat("*** Source file is set to           :", opt$source, "\n")
    cat("*** ROIs file is set to             :", opt$rois, "\n")    
    cat("*** Destination directory is set to :", opt$destination, "\n")
    cat("*** Window type is set to           :", opt$window,  "\n")
    cat("*** Window width is set to          :", opt$width,  "\n")
    cat("*** Window step is set to           :", opt$step,  "\n")
    if (opt$execute) 
        cat("*** Executing 3dNetCorr commands\n")    
    if ( ! is.null(opt$extra)) 
        cat("*** Extra agruments for 3dNetCorr   :", opt$extra,  "\n")
}


make.windows <- function (in.opt, in.length) {
    if (in.opt$window == "none" ) {
        windows=("0..$")
    } else if (in.opt$window == "consecutive" ) {

        window.starts=seq.int(from=0, to=(in.length-in.opt$width), by=in.opt$width)
        window.ends=(window.starts + in.opt$width) -1
        if (tail(window.ends, 1) < in.length) {
            if (! in.opt$quiet)
                cat("*** WARNING: The dataset is not a multiple of the window with and will be truncated\n")
        }
        
        indices=cbind(window.starts, window.ends)
        windows=unlist(apply(indices, 1, function (xx) {
            sprintf("%d..%d", xx[1], xx[2])
        }))

        return (windows)
        
        ## n.windows=ceiling(in.length/in.opt$width)
        ## windows=vector(mode="character", length=n.windows)
        ## cat("*** Preparing", n.windows, "windows: ")
        ## window.count=1
        ## done=FALSE
        ## window.start=0
        ## window.end=in.opt$width - 1
        ## while ( window.count <=  n.windows) {
        ##     if  (window.count ==  n.windows )
        ##         windows[window.count]=sprintf("%d..%d", window.start, in.length-1)
        ##     else 
        ##         windows[window.count]=sprintf("%d..%d", window.start, window.end)            
            
        ##     window.start=window.end + 1
        ##     window.end=window.end + in.opt$width

        ##     window.count=window.count+1
        ## }
        ## ## ww=seq.int(from=0, to=in.length, by=in.opt$width)
        ## cat(windows, sep=", ")
        ## cat("\n")
    } else { ## must be overlapping

        ## inspired by the following code
        ## http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r
        ## slideFunct <- function(data, window, step){
        ##     total <- length(data)
        ##     spots <- seq(from=1, to=(total-window), by=step)
        ##     result <- vector(length = length(spots))
        ##     for(i in 1:length(spots)){
        ##         result[i] <- mean(data[spots[i]:(spots[i]+window)])
        ##     }
        ##     return(result)
        ## }

        window.starts=seq.int(from=0, to=(in.length-in.opt$width), by=in.opt$step)
        window.ends=(window.starts + in.opt$width) -1
        if (tail(window.ends, 1) < in.length) {
            if (! in.opt$quiet)
                cat("*** WARNING: The dataset is not a multiple of the window with and will be truncated\n")
        }
        indices=cbind(window.starts, window.ends)
        windows=unlist(apply(indices, 1, function (xx) {
            sprintf("%d..%d", xx[1], xx[2])
        }))
    }
    return(windows)
}

make.3dNetCorr.commands <- function (in.opt, in.windows, in.sub.brik.selector) {

    yy=cbind(seq.int(1, length(in.windows)), in.windows)

    args=unlist(apply(yy, 1, function (xx) {
        if (in.sub.brik.selector) {
            sprintf("%s-prefix %s/%s.%02d -in_rois %s -inset %s\'[%s]\'",
                    ifelse(is.null(in.opt$extra), "", paste(in.opt$extra, " ", sep="")), in.opt$destination, in.opt$prefix, as.integer(xx[1]), in.opt$rois, in.opt$source, xx[2])
        } else {
            sprintf("%s-prefix %s/%s.%02d -in_rois %s -inset %s",
                    ifelse(is.null(in.opt$extra), "", paste(in.opt$extra, " ", sep="")), in.opt$destination, in.opt$prefix, as.integer(xx[1]), in.opt$rois, in.opt$source)
        }}))

    ## ll=unlist(sapply(in.windows, function (xx) {
    ##     sprintf("3dNetCorr -prefix %s/%s -rois %s -inset %s\'[%s]\'",
    ##             in.opt$destination, in.opt$prefix, in.opt$rois, in.opt$source, xx)
    ## }))
    
    ## ll=unlist(sapply(in.windows,
    ##     function (xx) {
    ##         sprintf("3dbucket -prefix %s.%s -session %s %s\'[%s]\'",
    ##                 in.opt$prefix, sub("..", ".to.", xx, fixed=TRUE), in.opt$destination, in.opt$source, xx)
    ##     }
    ##                  ))

    cmds=mapply(list, rep("3dNetCorr", length(args)), args, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    return(cmds)
    ## return(list("command" = "3dNetCorr", "arguments"=args))
}


read.mri.file <- function(in.filename) {
    if ( file.exists(in.filename)) {
        if (! opt$quiet)
            cat("*** Reading", in.filename, "\n")
        mri.dset=read.AFNI(in.filename, verb=TRUE)
    } else {
        mri.dset=NULL
        err.msg=paste("*** ERROR: No such file", in.filename, ". Cannot continue. Stopping.\n", sep="")
        if ( interactive()) {
            cat(err.msg)
        } else {
            stop(err.msg)
        }
    }
    return(mri.dset)
}


check.gridset <- function() { 
    ## check that the rois and the EPI are at the same resolution
    if ( ! all(abs(source.brik$delta) == abs(rois.brik$delta)) ) {
        err.msg="ERROR: The gridset of the source dataset and ROIs mask do not match. Cannot continue. Stopping.\n"
        if (interactive())
            cat(err.msg)
        else
            stop(err.msg)
    } else {
        return(TRUE)
    }
}

check.template.space <- function() {
    ## check that the rois and the EPI in the same space
    if (source.brik$NI_head$TEMPLATE_SPACE$dat != rois.brik$NI_head$TEMPLATE_SPACE$dat ) {
        err.msg=sprintf("ERROR: The template space of the source dataset (%s) and ROIs mask (%s) to not match. Cannot continue. Stopping.\n",
            source.brik$NI_head$TEMPLATE_SPACE$dat, rois.brik$NI_head$TEMPLATE_SPACE$dat)
        if (interactive())
            cat(err.msg)
        else
            stop(err.msg)
    } else {
        return(TRUE)
    }
}

conditionally.make.destination.dir <- function (in.opt) {
    if ( ! dir.exists(in.opt$destination)) {
        if (! in.opt$quiet)
            cat("*** Recursively creating", in.opt$destination, "\n")
        dir.create(in.opt$destination, recursive=TRUE)
    }
}


execute.commands <- function (in.commands) {
    for (ii in 1:length(in.commands)) {
        cat("*** Executing:", paste(net.corr.commands[[ii]][[1]],  net.corr.commands[[ii]][[2]]), "\n")
        ## the first element of each sublist in net.corr.commands
        ## contains the commands to be executed and the second the
        ## arguments to that command
        system2(net.corr.commands[[ii]][[1]], net.corr.commands[[ii]][[2]])
    }
}

print.commands <- function (in.commands) {
    for (ii in 1:length(in.commands)) {
        cat(paste(net.corr.commands[[ii]][[1]],  net.corr.commands[[ii]][[2]]), "\n")
    }
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

NO_ARGUMENT="0"
REQUIRED_ARGUMENT="1"
OPTIONAL_ARGUMENT="2"

## process command line arguments
spec = matrix(c(
    "help",          "h", NO_ARGUMENT,       "logical",   "Help for this program",
    "window",        "w", REQUIRED_ARGUMENT, "character", "Window type: consecutive, overlap, or none",
    "width",         "t", REQUIRED_ARGUMENT, "integer",   "Width of the window",
    "step",          "e", REQUIRED_ARGUMENT, "integer",   "Step between start of operlapping windows",    
    "source",        "s", NO_ARGUMENT,       "character", "Source file",
    "destination",   "d", REQUIRED_ARGUMENT, "character", "Destination directory",
    "rois",          "r", REQUIRED_ARGUMENT, "character", "Mask of ROIs",
    "prefix",        "p", REQUIRED_ARGUMENT, "character", "Prefix to use for created bucket files",
    "extra",         "x", REQUIRED_ARGUMENT, "character", "Extra arguments to provide to 3dNetCorr. These must be valid 3dNetCorr arguments",
    "quiet",         "q", NO_ARGUMENT,       "logical",   "Print no informational message. Only print 3dNetCorr commands.",
    "execute",       "u", NO_ARGUMENT,       "logical",   "Execute the 3dNetCorr commands instead of printing them out. The default is to print them out"
), byrow=TRUE, ncol=5)

if (interactive()) {
    cat("*** Setting interactive options\n")
    ## these are default arguments that are useful for testing
    ## purposes.
    ## args=c(
    ## "--source", "../data/425_A/rsfcPreprocessed/425_A.pm.cleanEPI.MNI.nii.gz",
    ## "-d", "../data/425_A/rsfcWindows",
    ## "-p", "425_A.pm.cleanEPI.window",
    ## "-r", "../standard/HarvardOxford/HarvardOxford-cort-maxprob-thr50-3mm.nii.gz",
    ## "-w", "consecutive",
    ## ## "-w", "none",
    ## "-w", "overlap",                
    ## "-t", "5",
    ## "-e", "3"#,
    ## "-x", "\"-part_corr\""
    ## )
    
    args=c(
        "--source", "../data/425_A/rsfcPreprocessed/425_A.pm.cleanEPI.MNI.nii.gz",
        "-d",  "../data/425_A/rsfcGraphs",
        "-p",  "425_A.pm.cleanEPI.aal2.whole.ts",
        "-w",  "none",
        "--extra", "\\-fish_z",
        "-r",  "../standard/aal2_for_SPM12/aal2.3mm.nii.gz",
        "-u"
    )
    ## args=c(
    ##     "--source", "../data/425_A/425_A.anat_struc_GM_to_template_GM_mod_s2.nii.gz",
    ##     "-d",  "../data/425_A/anatGraphs",
    ##     "-p",  "425_A.anat.aal2.whole",
    ##     "-w",  "none",
    ##     "-r",  "../standard/aal2_for_SPM12/aal2.nii.gz",
    ##     "-u"
    ## )
    opt = getopt(spec, opt=args)
} else {
    opt = getopt(spec)
}

opt=check.command.line.arguments(opt)
if (! opt$quiet)
    print.command.line.arguments.summary()

source.brik=read.mri.file(opt$source)
rois.brik=read.mri.file(opt$rois)

number.of.trs=source.brik$dim[4]
if (! opt$quiet)
    cat("*** EPI dataset is", number.of.trs, "TRs long.\n")

number.of.rois=max(max(rois.brik$brk))
if (! opt$quiet)
    cat("*** There are", number.of.rois, "ROIs in the ROI mask dataset.\n")


if(check.gridset())
    if (! opt$quiet)
        cat("*** Gridsets match\n")

if(check.template.space())
    if (! opt$quiet)
        cat("*** Template spaces match\n")

windows=make.windows(opt, number.of.trs)

conditionally.make.destination.dir(opt)

if (number.of.trs > 1) {
    net.corr.commands=make.3dNetCorr.commands(opt, windows, TRUE)
} else { 
    net.corr.commands=make.3dNetCorr.commands(opt, windows, FALSE)
}

if (opt$execute) {
    execute.commands(net.corr.commands)
} else {
    print.commands(net.corr.commands)
}
