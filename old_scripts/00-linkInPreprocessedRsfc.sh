#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM

## set -x 

root=/data/sanDiego/machLearnT1Analysis/data
rsfcRoot=/data/sanDiego/rsfcGraphAnalysis

subjects=$( cd $root; ls -1d [0-9][0-9][0-9]_A* )

#dirs=$( cd $root; ls -1d 117_A2BRIKS 334_A2BRIKS )

#dirs=$( cd $root; ls -1d 105_ABRIKS )

echo $dirs
subjectsDataDir=../data

for subject in $subjects ; do

    rsfcDir=${rsfcRoot}/data/${subject}/rsfcPreprocessed/
    
    if [[ -d ${rsfcDir} ]] ; then
	if [[ ! -h $root/$subject/rsfcPreprocessed ]]  ; then
	    echo "*** Linking in preprocessed RSFC data directory for subject ${subject}"
	    ( cd $root/$subject ; ln -sf ${rsfcDir}; symlinks -c ./ )
	else
	    echo "*** Skipping subject ${subject}: Link already exists"
	fi
    else
	echo "No such file ${rsfcDir}"
    fi
done
