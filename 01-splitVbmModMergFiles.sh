#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM

## set -x 

root=/data/sanDiego/machLearnT1Analysis/data

cd ../data/vbm
cd stats

mkdir split_mod_merg

cd split_mod_merg
## fslsplit ../GM_mod_merg_s2.nii.gz -t

cd ../../struc

files=$( ls *_GM_to_template_GM_mod.* )

## fileCount=$( ls *_GM_to_template_GM_mod.* | wc -l )

cd ../stats/split_mod_merg

modMergDir=$( pwd ) 
(( ii=0 ))
for modMergFile in ${files} ; do
    subjectNumber=$( echo $modMergFile | cut -d'.' -f 2 )
    vfile=$( printf "vol%04d.nii.gz" $ii )

    mm=$( echo $modMergFile | sed -e 's/....//' )
    targetDir=$root/$subjectNumber
    targetFile=${mm%%.nii.gz}_s2.nii.gz
    echo "*** Linking $modMergDir/$vfile to $targetDir/$targetFile"
    ln -sf $modMergDir/$vfile $targetDir/$targetFile
    ( cd $targetDir ; symlinks -c ./ )
    (( ii=ii+1 ))
done

