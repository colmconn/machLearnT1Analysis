#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

studyName=machLearnT1Analysis

GETOPT=$( which getopt )
ROOT=/data/sanDiego/$studyName
DATA=$ROOT/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

GETOPT_OPTIONS=$( $GETOPT  -o "s:l:p:r:" --longoptions "subject:,seedlist:,preproc:,rsfc:" -n ${programName} -- "$@" )
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
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	-p|--preproc)
	    preprocDir=$2; shift 2 ;;	    
	-r|--rsfc)
	    rsfcDir=$2; shift 2 ;;	    
	--) 
	    shift ; break ;;

	*) 
	    error_message "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $rsfcDir ] ; then 
    error_message "ERROR: The name of the direcotry into which to save the RSFC data was not provided. Exiting"
    exit
fi

if [ -z $preprocDir ] ; then 
    error_message "ERROR: The name of the direcotry where the preoprocessed RSFC data is stored was not provided. Exiting"
    exit
fi

if [ -z $subjectNumber ] ; then 
    error_message "ERROR: The subject ID was not provided. Exiting"
    exit
fi

if [ ! -f $seedList ] ; then
    error_message "ERROR: The seed list file does not exit. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList ) )
fi

info_message "Computing RSFC for the following seeds:"
info_message $seeds

if [[ ! -d $DATA/$subjectNumber/$rsfcDir ]] ; then
    mkdir $DATA/$subjectNumber/$rsfcDir
fi
cd $DATA/$subjectNumber/$rsfcDir

preprocessedRsfcDir=$DATA/$subjectNumber/$preprocDir

if [[ -f ${preprocessedRsfcDir}/errts.${subjectNumber}.tproject+tlrc.HEAD ]] ; then 
    for seed in $seeds ; do

	seedName=${seed##*/}
	if echo $seedName | grep -q "nii" ; then 
	    seedName=${seedName%%.nii*}
	else 
	    seedName=${seedName%%+*}
	fi

	if [[ ! -d ${seedName} ]] ; then
	    mkdir ${seedName}
	fi

	info_message "Extracting timeseries for seed ${seed}"
	3dROIstats -quiet -mask_f2short -mask ${seed} ${preprocessedRsfcDir}/errts.${subjectNumber}.tproject+tlrc.HEAD > ${seedName}/${seedName}.ts.1D

	info_message "Computing Correlation for seed ${seedName}"
	3dfim+ -input ${preprocessedRsfcDir}/errts.${subjectNumber}.tproject+tlrc.HEAD -ideal_file ${seedName}/${seedName}.ts.1D -out Correlation -bucket ${seedName}/${seedName}_corr
	
	info_message "Z-transforming correlations for seed ${seedName}"
	3dcalc -datum float -a ${seedName}/${seedName}_corr+tlrc.HEAD -expr 'log((a+1)/(a-1))/2' -prefix ${seedName}/${seedName}.z-score

	3drefit -sublabel 0 $subjectNumber ${seedName}/${seedName}.z-score+tlrc.HEAD
    done
fi
