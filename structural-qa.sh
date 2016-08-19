#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

DATA=/Volumes/data/sanDiego/machLearnT1Analysis/data

export AFNI_NOSPLASH=YES
#export AFNI_LAYOUT_FILE=elvis
#subjects="311_A"

## subjects="$*"


PIF=StructuralQa    #A string identifying programs launched by this script
                            #Get a free line and tag programs from this script
NPB="-npb `afni -available_npb_quiet` -pif $PIF -echo_edu" 

@Quiet_Talkers -pif $PIF > /dev/null 2>&1   #Quiet previously launched programs

# afni $NPB -niml -yesplugouts $adir/afni  >& /dev/null &

# plugout_drive  $NPB     

#for subject in $subjects ; do
while true ; do
    echo -n "Enter subject ID: "
    read subject
    
    if [[ "x$subject" == "x" ]] ; then
	break
    fi
    echo "####################################################################################################"
    echo "### Subject: $subject"
    echo "####################################################################################################"

    cd $DATA/$subject/
    afni $NPB -niml -yesplugouts  2> /dev/null &
    
    sleep 5
    #echo "Press enter for the original anatomy"
    #read
    plugout_drive $NPB \
	-com "SWITCH_UNDERLAY $subject.anat+orig.HEAD" \
	-com "SET_THRESHNEW 0" \
        -com 'OPEN_WINDOW A.axialimage geom=400x400+416+430 \
                     opacity=7'                         \
        -com 'OPEN_WINDOW A.coronalimage geom=+10+830 \
                     opacity=7'                         \
        -com 'OPEN_WINDOW A.sagittalimage geom=+10+430     \
                     opacity=7'                         \
	-quit

    echo "Enter clipping plane at the top of the brain"
    read above
    echo "Enter clipping plane at bottom of the brain"
    read bottom
    
    @clip_volume -input $subject.anat+orig.HEAD -above $above -below $bottom

    #echo "Press enter to quit afni and go to the next subject"
    #read
    plugout_drive $NPB \
	-com "QUIT" \
	-quit
    echo
    ## echo "Enter subject ID: "
done

echo "####################################################################################################"
echo "### All done!"
echo "####################################################################################################"
