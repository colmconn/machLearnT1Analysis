#!/bin/bash

## set -x 

cd /data/sanDiego/machLearnT1Analysis/data

subjects=$( ls -d [0-9][0-9][0-9]_[A]* )

## excessive motion threshold = 0.30
#LOG_TAIL=afniPreprocFromTemplate.log
#PREPROC_DIR=afniRsfcPreprocessed.noanaticor.NL.0.30

## excessive motion threshold = 0.25
LOG_TAIL=afniPreproc.noanaticor.NL.log
PREPROC_DIR=afniRsfcPreprocessed.noanaticor.NL.0.25

## subjects="105_A 106_A 107_A 108_A 109_A 111_A"

# subjects="303_A"

excessiveMotionThresholdFraction=0.2
echo "subject,anat HEAD,anat BRIK,resting HEAD,resting BRIK,all files OK?,resting nt,xmat?,numberOfCensoredVolumes,cutoff,inc/exc,sufficient data?"
for ss in $subjects ; do
    echo -n "$ss,"
    nnFiles=0
    if [[ -f $ss/${ss}.anat+orig.HEAD ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ -f $ss/${ss}.anat+orig.BRIK.gz ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi


    if [[ -f $ss/${ss}.resting+orig.HEAD ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ -f $ss/${ss}.resting+orig.BRIK.gz ]] ; then 
	echo -n "TRUE,"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "FALSE,"
    fi

    if [[ $nnFiles -eq 4 ]] ; then  
	echo -n "TRUE,"
    else 
	echo -n "FALSE,"
    fi
    
    if [[ -f $ss/${ss}.resting+orig.BRIK.gz ]] ; then 
	echo -n "$(3dinfo -nt $ss/${ss}.resting+orig.BRIK.gz ),"
	(( nnFiles=nnFiles + 1 ))
    else 
	echo -n "NA,"
    fi

    xmat_regress=X.xmat.1D 
    if [[ -f $ss/$PREPROC_DIR/$xmat_regress ]] ; then 
	fractionOfCensoredVolumes=$( 1d_tool.py -infile $ss/$PREPROC_DIR/$xmat_regress -show_tr_run_counts frac_cen )
	numberOfCensoredVolumes=$(   1d_tool.py -infile $ss/$PREPROC_DIR/$xmat_regress -show_tr_run_counts trs_cen )
	totalNumberOfVolumes=$(      1d_tool.py -infile $ss/$PREPROC_DIR/$xmat_regress -show_tr_run_counts trs_no_cen )

	## cutoff=$( echo "scale=0; $excessiveMotionThresholdFraction*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )
	## rounding method from http://www.alecjacobson.com/weblog/?p=256
	cutoff=$( echo "($(echo "scale=0;$excessiveMotionThresholdFraction*$totalNumberOfVolumes" | bc)+0.5)/1" | bc )
	echo -n "TRUE,$numberOfCensoredVolumes,$cutoff,"
	if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then
	    echo -n "EXCLUDE,"
	else
	    echo -n "INCLUDE,"	    
	fi
    else
	if [[ -f  $ss/${ss}-rsfc-${LOG_TAIL} ]] ; then 
	    censorCounts=( $( grep -a "Number of time points" $ss/${ss}-rsfc-${LOG_TAIL} | awk -F ";|:" '{print $2, $3}' | awk '{print $1, $4}') )
	    # ${censorCounts[0]} = total number of volumes
	    # ${censorCounts[1]} = total number of volumes left after censoring
	    totalNumberOfVolumes=${censorCounts[0]}
	    numberOfCensoredVolumes=$( expr ${censorCounts[0]} - ${censorCounts[1]} )
	    cutoff=$( echo "($(echo "scale=0;$excessiveMotionThresholdFraction*$totalNumberOfVolumes" | bc)+0.5)/1" | bc )
	    echo -n "FALSE,$numberOfCensoredVolumes,$cutoff,"
	    if [[ $numberOfCensoredVolumes -gt $cutoff ]] ; then
		echo -n "EXCLUDE,"
	    else
		echo -n "INCLUDE,"	    
	    fi
	else 
	    echo -n "FALSE,NA,NA,EXCLUDE,"
	fi
    fi

    if [[ -f  $ss/${ss}-rsfc-${LOG_TAIL} ]] ; then 
	if grep -qs "Insufficient data" $ss/${ss}-${LOG_TAIL} ; then 
	    echo "FALSE"
	else
	    echo "TRUE"
	fi
    else
	    echo "NLF"	
    fi

done


#{118_C,120_C,149_C,150_C,301_C,311_A,315_C,364_A,364_C}

#118_C 120_C 149_C 150_C 301_C 311_A 315_C 364_A 364_C
