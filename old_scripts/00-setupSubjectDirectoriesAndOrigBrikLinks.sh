#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM
set -x 

root=/data/sanDiego/
dirs=$( cd $root; ls -1d [0-9][0-9][0-9]_ABRIKS )
#dirs=$( cd $root; ls -1d 117_A2BRIKS 334_A2BRIKS )

#dirs=$( cd $root; ls -1d 105_ABRIKS )

echo $dirs
subjectsDataDir=../data

for dir in $dirs ; do

    subject=${dir%%BRIKS}
    if [[ -f $root/$dir/$subject+orig.HEAD ]] ; then
	(cd $subjectsDataDir; mkdir $subject )
	(cd $subjectsDataDir/$subject; ln -sf $root/$dir/${subject}+orig.HEAD . )
	(cd $subjectsDataDir/$subject; ln -sf $root/$dir/${subject}+orig.BRIK.gz . )
	(cd $subjectsDataDir/$subject; symlinks -c ./ )
	

    else
	echo "No such file $root/$dir/$subject+orig.HEAD"
    fi

done
