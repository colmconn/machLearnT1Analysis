#!/bin/bash

## set -x 

programName=`basename $0`

trap exit SIGHUP SIGINT SIGTERM

studyName=machLearnT1Analysis

GETOPT=$( which getopt )
ROOT=/data/sanDiego/$studyName
DATA=$ROOT/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

function doZeropad {
    local subject="$1"
    if [[ $subject == "341_A" ]] ; then
	sup="-S 30"
    fi
    info_message "Zeropadding anat and EPI for subject $subject"
    if [[ -f $DATA/$subject/${subject}.anat_clp+orig.HEAD ]] ; then
	if [[ $force -eq 1 ]] || \
	   [[ ! -f $DATA/$subject/${subject}.anat.zp+orig.HEAD ]]  || \
	   [[ $DATA/$subject/${subject}.anat_clp+orig.HEAD -nt $DATA/$subject/${subject}.anat.zp+orig.HEAD ]] ; then
	    ( cd $DATA/$subject ; 3dZeropad -I 30 $sup -prefix ${subject}.anat.zp ${subject}.anat_clp+orig.HEAD )
	fi
    else
	if [[ $force -eq 1 ]] || \
	   [[ ! -f $DATA/$subject/${subject}.anat.zp+orig.HEAD ]] || \
	   [[ $DATA/$subject/${subject}.anat+orig.HEAD -nt $DATA/$subject/${subject}.anat.zp+orig.HEAD ]]; then 
	    ( cd $DATA/$subject ; 3dZeropad -I 30 $sup -prefix ${subject}.anat.zp ${subject}.anat+orig.HEAD )
	fi
    fi
    if [[ $force -eq 1 ]] || [[ ! -f $DATA/$subject/${subject}.resting.zp+orig.HEAD ]] ; then 
	( cd $DATA/$subject ; 3dZeropad -I 30 $sup -prefix ${subject}.resting.zp ${subject}.resting+orig.HEAD )
    fi
}

GETOPT_OPTIONS=$( $GETOPT \
		      -o "fe:m:o:h:l:h:b:t:nq" \
		      --longoptions "force,excessiveMotionThresholdFraction:,motionThreshold:,outlierThreshold:,threads:,lowpass:,highpass:,blur:,tcat:,nonlinear,enqueue" \
		      -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

## 1 = force creation of zero padded files
force=0

## enqueue the job for execution
enqueue=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-f|--force)
	    force=1; shift 1;;
	-e|--excessiveMotionThresholdFraction)
	    excessiveMotionThresholdFraction=$2; shift 2 ;;	
	-m|--motionThreshold)
	    motionThreshold=$2; shift 2 ;;	
	-o|--outlierThreshold)
	    outlierThreshold=$2; shift 2 ;;	
	-h|--threads)
	    threads=$2; shift 2 ;;	
	-l|--lp)
	    lowpass=$2; shift 2 ;;	
	-h|--hp)
	    highpass=$2; shift 2 ;;	
	-b|--blur)
	    blur=$2; shift 2 ;;	
	-t|--tcat)
	    tcat=$2; shift 2 ;;	
	-n|--nonlinear)
	    nonlinear=1; shift 1 ;;	
	-q|--enqueue)
	    enqueue=1; shift 1 ;;	
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [[ $force -eq 1 ]] ; then
    info_message "Forcing recreation of ZEROPADed files"
fi

####################################################################################################
## Check that appropriate values are used to initialize arguments that
## control analysis if no values were provided on the command line

## The following values are used to exclude subjects based on the
## number of volumes censored during analysis
if [[ "x$excessiveMotionThresholdFraction" == "x" ]] ; then
    excessiveMotionThresholdFraction=0.2
    excessiveMotionThresholdPercentage=20
    warn_message "No excessiveMotionThresholdFraction threshold was provided. Defaulting to $excessiveMotionThresholdFraction => ${excessiveMotionThresholdPercentage}%"
else
    excessiveMotionThresholdPercentage=$( echo "(($excessiveMotionThresholdFraction*100)+0.5)/1" | bc ) 

    info_message "Using ${excessiveMotionThresholdFraction} as the subject exclusion motion cutoff fraction"
    info_message "Using ${excessiveMotionThresholdPercentage}% as subject exclusion motion cutoff percentage"
    info_message "Note that these values are used to exclude subjects based on the number of volumes censored during analysis"
fi


## motionThreshold and outlierThreshold are the values passed to
## afni_proc.py and are used when deciding to censor a volume or not
if [[ "x${motionThreshold}" == "x" ]] ; then
    motionThreshold=0.2
    warn_message "No motionThreshold value was provided. Defaulting to $motionThreshold"
else
    info_message "Using motionThreshold of ${motionThreshold}"
fi

if [[ "x${outlierThreshold}" == "x" ]] ; then
     outlierThreshold=0.1
     warn_message "No outlierThreshold value was provided. Defaulting to $outlierThreshold"
else
    info_message "Using outlierThreshold of ${outlierThreshold}"
fi

if [[ "x${threads}" == "x" ]] ; then
     threads=1
     warn_message "No value for the number of parallel threads to use was provided. Defaulting to $threads"
else
    info_message "Using threads value of ${threads}"
fi

if [[ "x${lowpass}" == "x" ]] ; then
     lowpass="0.01"
     warn_message "No value for lowpass filter value to use was provided. Defaulting to $lowpass"
else
    info_message "Using lowpass filter value of ${lowpass}"
fi

if [[ "x${highpass}" == "x" ]] ; then
     highpass="0.1"
     warn_message "No value for highpass filter value to use was provided. Defaulting to $highpass"
else
    info_message "Using highpass filter value of ${highpass}"
fi

if [[ "x${blur}" == "x" ]] ; then
     blur="8"
     warn_message "No value for blur filter value to use was provided. Defaulting to $blur"
else
    info_message "Using blur filter value of ${blur}"
fi

if [[ "x${tcat}" == "x" ]] ; then
     tcat="3"
     warn_message "No value for tcat, the number of TRs to censor from the start of each volume, was provided. Defaulting to $tcat"
else
    info_message "Using tcat filter value of ${tcat}"
fi

if [[ $nonlinear -eq 1 ]] ; then 
    info_message "Using nonlinear alignment"
    scriptExt="noanaticor.NL"
else 
    info_message "Using affine alignment only"
    scriptExt="noanaticor.aff"    
fi

####################################################################################################
if [[ "$#" -gt 0 ]] ; then
    subjects="$@"
else
    subjects=$( cd $DATA ; ls -d *_A* )
fi

## subjects="107_A 109_A 137_A 151_A 153_A 154_A 164_A 300_A 304_A 325_A 332_A 342_A 343_A 345_A 348_A 351_A 356_A 360_A 366_A 378_A 376_A 380_A 395_A"

## subjects="105_A"

[[ -d run ]] || mkdir run

for subject in $subjects ; do
    info_message "#################################################################################################"
    info_message "Generating script for subject $subject"

    if  [[ ! -f $DATA/$subject/${subject}.resting+orig.HEAD ]] && \
	[[ ! -f $DATA/$subject/${subject}.resting+orig.BRIK.gz ]]  ; then

	warn_message "Can not find resting state EPI file for ${subject}. Skipping."
	continue
    else
	epiFile=$DATA/$subject/${subject}.resting+orig.HEAD
    fi

    if  [[ ! -f $DATA/$subject/${subject}.anat+orig.HEAD ]] && \
	[[ ! -f $DATA/$subject/${subject}.anat+orig.BRIK.gz ]]  ; then

	warn_message "Can not find anatomy file for subject ${subject}. Skipping."
	continue
    else
	anatFile=${DATA}/$subject/$subject.anat+orig.HEAD
    fi

    if [[ $nonlinear -eq 1 ]] ; then 
	outputScriptName=run/run-afniRsfcPreproc-${subject}.${scriptExt}.sh
    else
	outputScriptName=run/run-afniRsfcPreproc-${subject}.${scriptExt}.sh	
    fi
    
    case $subject in
	####################################################################################################
	## first batch of alignment assessments
    	109_A)
    	    extraAlignmentArgs="-align_opts_aea  -cost lpa -giant_move"
    	    ;;
    	113_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	120_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	133_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	135_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	137_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc -giant_move"
    	    ;;
    	143_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	145_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	153_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpa -giant_move"
    	    ;;
    	154_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	155_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc -giant_move"
    	    ;;
    	159_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc -giant_move"
    	    ;;
    	165_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpa -giant_move"
    	    ;;
    	307_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpa -giant_move"
    	    ;;
    	308_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	311_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	313_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc -giant_move"
    	    ;;
    	320_A)
    	    ## use_clipped_anat
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc -giant_move"
    	    ;;
    	325_A)
    	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea  -cost lpc+ZZ -giant_move"
    	    ;;
    	330_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	334_A)
   	    doZeropad $subject	    
   	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
   	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
    	341_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	345_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;	
    	346_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	347_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	348_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
    	356_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;	
    	364_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
    	370_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	374_A)
    	    # doZeropad $subject	    
    	    # anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    # epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    # extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	380_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
    	390_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	411_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	414_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;	
    	411_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
	####################################################################################################
	## second batch of alignment assessments
	108_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
	118_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
	116_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	121_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	124_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	126_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	127_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	132_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	142_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	146_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	150_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	160_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	164_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	167_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	300_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;	
    	306_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;	
    	310_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	312_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	314_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	315_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	317_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	321_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	322_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	326_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	331_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	333_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	335_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	338_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	343_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	351_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	357_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	360_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	362_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	367_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	374_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	386_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	389_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	391_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	395_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	397_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	398_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	403_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;
    	406_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
    	415_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	415_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	417_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	421_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	423_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	424_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
	####################################################################################################
	## 3rd batch
    	339_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	369_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	371_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	372_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	373_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	377_A)
    	    doZeropad $subject	    
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	376_A)
    	    anatFile=${DATA}/$subject/$subject.anat+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
	####################################################################################################
	## 4th batch
    	117_A2)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
	####################################################################################################
	## 5th batch
    	106_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	114_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	119_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ"
    	    ;;
    	132_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	158_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	165_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;
    	167_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	307_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc -giant_move"
    	    ;;
    	316_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
    	386_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ"
    	    ;;	    
    	424_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc"
    	    ;;
	####################################################################################################
	## 6th batch
    	300_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpa -giant_move"
    	    ;;	
    	336_A)
	    doZeropad $subject
    	    anatFile=${DATA}/$subject/$subject.anat.zp+orig.HEAD
    	    epiFile=$DATA/$subject/${subject}.resting.zp+orig.HEAD
    	    extraAlignmentArgs="-align_opts_aea -cost lpc+ZZ -giant_move"
    	    ;;	
	*)
    	    extraAlignmentArgs=""
    	    ;;
    esac

    ## do non-linear warping? If so add the flag to the extra
    ## alignment args variable
    if [[ $nonlinear -eq 1 ]] ; then 
	extraAlignmentArgs="${extraAlignmentArgs} -tlrc_NL_warp"
	anat_base=$( basename $anatFile )
	anat_base=${anat_base%%+*}
	if [[ -f ${DATA}/${subject}/afniRsfcPreprocessed.NL/${anat_base}_al_keep+tlrc.HEAD ]] && \
	   [[ -f ${DATA}/${subject}/afniRsfcPreprocessed.NL/anat.un.aff.Xat.1D ]] && \
	   [[ -f ${DATA}/${subject}/afniRsfcPreprocessed.NL/anat.un.aff.qw_WARP.nii ]] ; then
	    info_message "Supplying prexisting nonlinear warped anatomy to afni_proc.py"
	    extraAlignmentArgs="${extraAlignmentArgs} \\
	     -tlrc_NL_warped_dsets ${DATA}/${subject}/afniRsfcPreprocessed.NL/${anat_base}_al_keep+tlrc.HEAD \\
                                   ${DATA}/${subject}/afniRsfcPreprocessed.NL/anat.un.aff.Xat.1D \\
                                   ${DATA}/${subject}/afniRsfcPreprocessed.NL/anat.un.aff.qw_WARP.nii"
	fi
    fi

    info_message "Writing script: $outputScriptName"


    cat <<EOF > $outputScriptName
#!/bin/bash

set -x 

#$ -S /bin/bash

## disable compression of BRIKs/nii files
unset AFNI_COMPRESSOR

export PYTHONPATH=$AFNI_R_DIR

## use the newer faster despiking method. comment this out to get the
## old one back
export AFNI_3dDespike_NEW=YES

# turn off anoying colorization of info/warn/error messages since they
# only result in gobbledygook
export AFNI_MESSAGE_COLORIZE=NO

## only use a single thread since we're going to run so many subjects
## in parallel
export OMP_NUM_THREADS=${threads}

excessiveMotionThresholdFraction=$excessiveMotionThresholdFraction
excessiveMotionThresholdPercentage=$excessiveMotionThresholdPercentage

cd $DATA/$subject

preprocessingScript=${subject}.afniRsfcPreprocess.$scriptExt.csh
rm -f \${preprocessingScript}

outputDir=afniRsfcPreprocessed.$scriptExt
rm -fr \${outputDir}

motionThreshold=${motionThreshold}
outlierThreshold=${outlierThreshold}

##	     -tcat_remove_first_trs ${tcat}					\\
## -tlrc_opts_at -init_xform AUTO_CENTER \\
## 	     -regress_censor_outliers \$outlierThreshold                 	\\

afni_proc.py -subj_id ${subject}						\\
             -script \${preprocessingScript}					\\
	     -out_dir \${outputDir}						\\
	     -blocks despike tshift align tlrc volreg mask blur regress		\\
	     -copy_anat $anatFile                                               \\
	     -dsets $epiFile                                                    \\
	     -tlrc_base MNI_caez_N27+tlrc					\\
	     -volreg_align_to first    						\\
	     -volreg_tlrc_warp ${extraAlignmentArgs}				\\
	     -blur_size ${blur}                                                 \\
	     -blur_to_fwhm  							\\
	     -blur_opts_B2FW "-ACF -rate 0.2 -temper"                           \\
	     -mask_apply group							\\
	     -mask_segment_anat yes						\\
	     -regress_censor_first_trs ${tcat}					\\
	     -mask_segment_erode yes						\\
	     -regress_ROI WMe							\\
	     -regress_bandpass ${lowpass} ${highpass}				\\
	     -regress_apply_mot_types demean   					\\
             -regress_censor_motion \$motionThreshold              		\\
	     -regress_run_clustsim no						\\
	     -regress_est_blur_epits                                            \\
	     -regress_est_blur_errts

if [[ -f \${preprocessingScript} ]] ; then 
   tcsh -xef \${preprocessingScript}

    cd \${outputDir}
    xmat_regress=X.xmat.1D 

    if [[ -f \$xmat_regress ]] ; then 

        fractionOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts frac_cen )
        numberOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_cen )
        totalNumberOfVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_no_cen )

        ## rounding method from http://www.alecjacobson.com/weblog/?p=256
        cutoff=\$( echo "((\$excessiveMotionThresholdFraction*\$totalNumberOfVolumes)+0.5)/1" | bc )
	if [[ \$numberOfCensoredVolumes -gt \$cutoff ]] ; then 

	    echo "*** A total of \$numberOfCensoredVolumes of
	    \$totalNumberOfVolumes volumes were censored which is
	    greater than \$excessiveMotionThresholdFraction
	    (n=\$cutoff) of all total volumes of this subject" > \\
		00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt

	    echo "*** WARNING: $subject will not be analysed due to having more than \${excessiveMotionThresholdPercentage}% of their volumes censored."
	fi

    	# trs=\$( 1d_tool.py -infile \$xmat_regress -show_trs_uncensored encoded   \\
        #                   -show_trs_run 01 )
    	# if [[ \$trs != "" ]] ; then  
	#     3dFWHMx -ACF -detrend -mask mask_group+tlrc                      \\
        # 	errts.$subject.anaticor+tlrc"[\$trs]" > acf.blur.errts.1D
	# fi
    else
	touch 00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt
    fi
    echo "Compressing BRIKs and nii files"
    find ./ \( -name "*.BRIK" -o -name "*.nii" \) -print0 | xargs -0 gzip
else
    echo "*** No such file \${preprocessingScript}"
    echo "*** Cannot continue"
    exit 1
fi	

EOF

    chmod +x $outputScriptName
    if [[ $enqueue -eq 1 ]] ; then
	info_message "Submitting job for execution to queuing system"
	LOG_FILE=$DATA/$subject/$subject-rsfc-afniPreproc.${scriptExt}.log
	info_message "To see progress run: tail -f $LOG_FILE"
	rm -f ${LOG_FILE}
	qsub -N rsfc-$subject -q all.q -j y -m n -V -wd $( pwd )  -o ${LOG_FILE} $outputScriptName
    else
	info_message "Job *NOT* submitted for execution to queuing system"
	info_message "Pass -q or --enqueue options to this script to do so"	
    fi

done

if [[ $enqueue -eq 1 ]] ; then 
    qstat
fi
