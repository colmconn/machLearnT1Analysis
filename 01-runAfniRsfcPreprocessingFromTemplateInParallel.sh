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

GETOPT_OPTIONS=$( $GETOPT \
		      -o "e:h:" \
		      --longoptions "excessiveMotionThresholdFraction:threads:" \
		      -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-e|--excessiveMotionThresholdFraction)
	    excessiveMotionThresholdFraction=$2; shift 2 ;;	
	-h|--threads)
	    threads=$2; shift 2 ;;	
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

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

if [[ "x${threads}" == "x" ]] ; then
     threads=1
     warn_message "No value for the number of parallel threads to use was provided. Defaulting to $threads"
else
    info_message "Using threads value of ${threads}"
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

    outputScriptName=run/run-afniRsfcPreprocFromTemplate-${subject}.sh

    info_message "Writing script: $outputScriptName"
    
    cat <<EOF > $outputScriptName
#!/bin/bash

set -x 

#$ -S /bin/bash

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

preprocessingScript=$SCRIPTS_DIR/afniRsfcPreprocessTemplate.csh

outputDir=afniRsfcPreprocessed
if [[ -d \$outputDir ]] && [[ ! -d \${outputDir}.orig ]] ; then
    cp -ra \$outputDir \${outputDir}.orig
fi

if [[ -f \${preprocessingScript} ]] ; then 
   tcsh -xef \${preprocessingScript} ${subject}

    cd \${outputDir}
    xmat_regress=X.xmat.1D
    ## always delete any preexisting file so we can start fresh each time this is executed
    rm -f 00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt
    if [[ -f \$xmat_regress ]] ; then 

        fractionOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts frac_cen )
        numberOfCensoredVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_cen )
        totalNumberOfVolumes=\$( 1d_tool.py -infile \$xmat_regress -show_tr_run_counts trs_no_cen )

	## cutoff=\$( echo "scale=0; \$excessiveMotionThresholdFraction*\$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )
        ## rounding method from http://www.alecjacobson.com/weblog/?p=256
        cutoff=\$( echo "(\$(echo "scale=0;\$excessiveMotionThresholdFraction*\$totalNumberOfVolumes" | bc)+0.5)/1" | bc )
	if [[ \$numberOfCensoredVolumes -gt \$cutoff ]] ; then 

	    echo "*** A total of \$numberOfCensoredVolumes of
	    \$totalNumberOfVolumes volumes were censored which is
	    greater than \$excessiveMotionThresholdFraction
	    (n=\$cutoff) of all total volumes of this subject" > \\
		00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt

	    echo "*** WARNING: $subject will not be analysed due to having more than \${excessiveMotionThresholdPercentage}% of their volumes censored."
	fi

    	trs=\$( 1d_tool.py -infile \$xmat_regress -show_trs_uncensored encoded   \\
                          -show_trs_run 01 )
    	if [[ \$trs != "" ]] ; then  
	    3dFWHMx -ACF -detrend -mask mask_group+tlrc                      \\
        	errts.$subject.anaticor+tlrc"[\$trs]" > acf.blur.errts.1D
	fi
    else
	touch 00_DO_NOT_ANALYSE_${subject}_\${excessiveMotionThresholdPercentage}percent.txt
    fi
else
    echo "*** No such file \${preprocessingScript}"
    echo "*** Cannot continue"
    exit 1
fi	

EOF

    chmod +x $outputScriptName
    LOG_FILE=$DATA/$subject/$subject-rsfc-afniPreprocFromTemplate.log
    rm -f ${LOG_FILE}
    qsub -N rsfc-$subject -q all.q -j y -m n -V -wd $( pwd )  -o ${LOG_FILE} $outputScriptName

done

qstat
