#!/bin/bash

ROOT=${MDD_ROOT:-/data/sanDiego/machLearnT1Analysis}
DATA=$ROOT/data

function makeMniLink {
    local directory=$1

    echo "*** Making MNI links in $directory"
    
    owd=$( pwd ) 
    cd $directory

    if [[ ! -f MNI_caez_N27+tlrc.HEAD ]] ; then
	for ff in ../../MNI_caez* ; do
	    ln -sf $ff
	done
    fi

    cd $owd
}

groupResultsDirs="$( ls -d $DATA/group.results/0.*/ttest.* )"

for dd in ${groupResultsDirs} ; do
    makeMniLink $dd
done

    
