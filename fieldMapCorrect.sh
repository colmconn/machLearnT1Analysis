#!/bin/bash

# set -x 

trap exit SIGHUP SIGINT SIGTERM

studyName=machLearnT1Analysis

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/$studyName}

DATA=$ROOT/data
RAW_DATA=/data/sanDiego
PROCESSED_DATA=$DATA
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

## ensure that ppge4 can find epidewarp.ucsd.sh
export PATH=$PATH:$SCRIPTS_DIR

GETOPT_OPTIONS=$( $GETOPT  -o "fs:" --longoptions "force,subject::" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    error_message "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

force=0
# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-f|--force)
	    force=1; shift ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    error_message "The subject ID was not provided."
    exit
fi

function do_fieldmap_correction {
    subject=$1
    dicomTask="$2"
    task="$3"
    
    subjectDicomContainerDir=$RAW_DATA/$subject
    if [[ ${#subjectDicomContainerDir} -gt 0 ]] && [[ -d $subjectDicomContainerDir ]] ; then 
	if [[ ! -f $subjectDicomContainerDir/sinfo.txt ]] ; then 
	    (cd $subjectDicomContainerDir; ../scheck )
	fi

	if [[ -d $subjectDicomContainerDir/asl ]] ; then
	    if [[ ! -f $subjectDicomContainerDir/asl/sinfo.txt ]] ; then
		(cd $subjectDicomContainerDir/asl; ../../scheck )
	    fi
	fi

	epiDirs=( $( cat $subjectDicomContainerDir/sinfo.txt $subjectDicomContainerDir/asl/sinfo.txt 2> /dev/null | grep -i "$dicomTask" | uniq | tail -1 | awk '{print $1}' ) )
	if [[ ${#epiDirs[@]} -gt 1 ]] ; then
	    error_message "This script cannot currently handle fieldmap correction for moe than one run"
	    exit 1
	fi
	
	fmTe1Dir=$( cat $subjectDicomContainerDir/sinfo.txt $subjectDicomContainerDir/asl/sinfo.txt 2> /dev/null | grep -v DTI | grep -i "fm_TE1_NFS" | uniq | tail -1 | awk '{print $1}' )
	fmTe2Dir=$( cat $subjectDicomContainerDir/sinfo.txt $subjectDicomContainerDir/asl/sinfo.txt 2> /dev/null | grep -v DTI | grep -i "fm_TE2_NFS" | uniq | tail -1 | awk '{print $1}' )

	if [[ "x$fmTe1Dir" != "x" ]] && [[ "x$fmTe2Dir" != "x" ]] ; then
	    session=$PROCESSED_DATA/${subject}

	    if [[ -f $session/${subject}.${task}.nii.gz ]] ; then
		n_arg="-n $session/${subject}.${task}.nii.gz"
	    else
		n_arg=""
	    fi
	    ## fmc=field map corrected
	    ( cd $subjectDicomContainerDir;  $SCRIPTS_DIR/ppge4.sh  -d1 $fmTe1Dir -d2 $fmTe2Dir -i ${epiDirs[0]} -o "$session/${subject}.${task}.fmc" $n_arg )
	else
	    if [[ "x$fmTe1Dir" == "x" ]]; then 
		warn_message "No s-directories found for the field map TE1 in the sinfo.txt file. Cannot fieldmap correct $task task data for subject ${subject}. Skipping."
	    fi
	    if [[ "x$fmTe2Dir" == "x" ]]; then 
		warn_message "No s-directories found for the field map TE2 in the sinfo.txt file. Cannot fieldmap correct $task task data for subject ${subject}. Skipping."
	    fi
	fi
    else 
	warn_message "Cannot find $subjectDicomContainerDir"
	warn_message "Skipping"
    fi
}

##i=$((${#subjectNumber}-1))
##timepoint="${subjectNumber:$i:1}"

info_message "####################################################################################################"
info_message "### Subject: $subjectNumber"


if [[ ! -d $PROCESSED_DATA/${subjectNumber} ]] ; then
    info_message "Making the subject's directory in the study data folder: $PROCESSED_DATA/${subjectNumber}"
    mkdir $PROCESSED_DATA/${subjectNumber}
fi

dicomTask="resting state"
task="resting"

if [[ -f $PROCESSED_DATA/${subjectNumber}/${subjectNumber}.${task}+orig.HEAD ]] && [[ ! -f $PROCESSED_DATA/${subjectNumber}/${subjectNumber}.${task}.nii.gz ]] ; then
    ( cd $PROCESSED_DATA/${subjectNumber}; 3dcopy ${subjectNumber}.${task}+orig.HEAD ${subjectNumber}.${task}.nii )
fi
#if [[ $force -eq 1 ]] || [[ ! -f $PROCESSED_DATA/${subjectNumber}/${subjectNumber}.${task}+orig.HEAD ]] ; then 
    info_message "****************************************************************************************************"
    info_message "*** Resting state reconstruction"
    export AFNIDIR=${AFNI_R_DIR} 
    do_fieldmap_correction "$subjectNumber"  "$dicomTask" "$task"
#fi
