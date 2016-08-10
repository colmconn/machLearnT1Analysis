#!/bin/bash

## set -x 

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

GETOPT_OPTIONS=$( $GETOPT  -o "s:" --longoptions "subject::" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    error_message "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
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

function runDimon {
    dicomTask="$1"
    task="$2"
    run="$3"
    subjectDicomContainerDir="$4"
    sdir="$5"
    output="$6"
    
    ## dcmFile=${dcmfiles[$ii]}
    session=$PROCESSED_DATA/${subject}
    
    if [[ ${#run} -gt 0 ]] ; then
        ntask="$task$run"
    else
        ntask="$task"		    
    fi

    prefix="$subject.$ntask"
    info_message "Now creating AFNI HEAD/BRIK of the ${dicomTask} task for $subject"
    info_message "Prefix will be $prefix"
    if [[ $output == "afni" ]] || [[ $output == "both" ]] ; then 
	( cd $subjectDicomContainerDir;  \
	  Dimon -infile_pattern $sdir/'i*' -dicom_org -gert_filename make.$ntask \
		-save_file_list "${subject}.${ntask}.dicom.list" \
		-gert_to3d_prefix "${subject}.${ntask}" -gert_outdir ${session} \
			   -GERT_Reco  -gert_create_dataset -quit )
    fi
    # if [[ $output == "mgz" ]] || [[ $output == "both" ]] ; then   
    #     ( mkdir -p $session/../mri/orig/ ; 
    #       cd $subjectDicomContainerDir;  \
    #       mri_convert -it dicom -ot mgz -i $dcmFile -o $session/../mri/orig/001.mgz ) 
    # fi 
}


function reconstruct {
    output="$1"
    subject=$2
    dicomTask="$3"
    task="$4"
    
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
	
	sdirs=( $( cat $subjectDicomContainerDir/sinfo.txt $subjectDicomContainerDir/asl/sinfo.txt 2> /dev/null | grep -i "$dicomTask" | uniq | tail -1 | awk '{print $1}' ) )
	## dcmfiles=( $( cat $subjectDicomContainerDir/sinfo.txt | grep -i "$dicomTask" | tail -1 | awk '{print $2}' ) )

	if [[ ${#sdirs[@]} -gt 0 ]] ; then
	    (( run=1 ))
	    (( ii=0 ))
	    ## csdir = candidate sdir
	    for csdir in ${sdirs[$ii]} ; do
		if [[ ${#sdirs[@]} -eq 1 ]] ; then
		    ## if there is only 1 run set this to the empty string
		    ## so that the output file from Dimon is not tagged
		    ## with a useless numeric value
		    run=""
		fi
		if [[ -d  $subjectDicomContainerDir/$csdir ]] ; then
		    runDimon "$dicomTask" "$task" "$run" $subjectDicomContainerDir $csdir "afni"
		elif [[ -d  $subjectDicomContainerDir/asl/$csdir ]] ; then
		    runDimon "$dicomTask" "$task" "$run" $subjectDicomContainerDir/asl $csdir "afni"
		else
		    warn_message "Neither $subjectDicomContainerDir/$sdir nor $subjectDicomContainerDir/asl/$sdir worked for subject $subject. Fix this manually"
		fi
		(( run=run+1 ))
		(( ii=ii+1 ))
	    done
	else
	    warn_message "No s-directories found found in the sinfo.txt file. Cannot reconstruct $task task data for subject ${subject}. Skipping."
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


info_message "****************************************************************************************************"
info_message "T1 Anatomy reconstruction"

dicomTask="FSPGR_SAG_TI550"
task="anat"
if [[ ! -f $PROCESSED_DATA/${subjectNumber}/${subjectNumber}.${task}+orig.HEAD ]] ; then 
    ## reconstruct "mgz" "$subjectNumber"  "$dicomTask" "$task"
    reconstruct "afni" "$subjectNumber"  "$dicomTask" "$task"     
fi


info_message "****************************************************************************************************"
info_message "*** Resting state reconstruction"

dicomTask="resting state"
task="resting"
if [[ ! -f $PROCESSED_DATA/${subjectNumber}/${subjectNumber}.${task}+orig.HEAD ]] ; then 
    reconstruct "afni" "$subjectNumber"  "$dicomTask" "$task"
fi
