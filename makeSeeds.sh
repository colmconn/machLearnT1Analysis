#!/bin/bash

set -x 

trap exit SIGHUP SIGINT SIGTERM


[[ ! -d ../data/seeds ]] && mkdir ../data/seeds

cd ../data/seeds


## get the MNI_caez_N27 brain and convert to nii format
3dcopy /data/software/afni/MNI_caez_N27+tlrc.HEAD ./MNI_caez_N27.nii

mkdir MNI152_space

cd MNI152_space

## get the sgACC seeds
cp /data/sanDiego/rsfcGraphAnalysis/data/seeds/ACC/acci{8,9}L.nii.gz .
cp /data/sanDiego/rsfcGraphAnalysis/data/seeds/ACC/acci{8,9}R.nii.gz .


## get the amygdala and anterior inssula seeds 
for ff in /data/sanDiego/r01_preliminary/data/seeds/a*tlrc.HEAD ; do
    dd=${ff##*/}
    3dcopy $ff ./${dd%%+*}.nii
done


## make sure they overlay the MNI_caez_N27 (this step is probably not
## necessary since everything is in MNI space already

# for ff in a* ; do
#     flirt -in $ff -ref MNI_caez_N27.nii.gz -usesqform -applyxfm -out ${ff##*/} -interp nearestneighbour
# done

## ensure that all of the ROIs are on the same gridset as the the EPI
## timeseries. That is make them match the /MNI_caez_N27 brain
for ff in a*.nii.gz ; do
    3dresample -dxyz 3 3 3 -master ../MNI_caez_N27.nii.gz -inset $ff -prefix CAEZ.${ff%%.gz}
    mv -f CAEZ.${ff} ../$ff
done


echo "-42 16 36" | 3dUndump -master ./MNI_caez_N27.nii.gz -srad 8 -prefix L_DLPFC -xyz -
echo " 42 16 36" | 3dUndump -master ./MNI_caez_N27.nii.gz -srad 8 -prefix R_DLPFC -xyz -
for ff in *_DLPFC+tlrc.HEAD ; do
    3dresample -dxyz 3 3 3 -master ./MNI_caez_N27.nii.gz  -inset $ff -prefix ${ff%%+*}.3mm
done


# ## not copy the DLPFC and posterior cingulate seeds as Felipe has already made them
# cd ../
# cp /data/jain/seeds/L_DLPFC/L_DLPFC_3mm+tlrc.HEAD /data/jain/seeds/L_DLPFC/L_DLPFC_3mm+tlrc.BRIK ./
# cp /data/jain/seeds/R_DLPFC/R_DLPFC_3mm+tlrc.HEAD /data/jain/seeds/R_DLPFC/R_DLPFC_3mm+tlrc.BRIK ./

# ## now create a bilateral version
# 3dcalc -a L_DLPFC_3mm+tlrc.HEAD -b R_DLPFC_3mm+tlrc.HEAD -expr "step(a) + step(b)" -prefix DLPFC.bilateral.nii

# cp /data/jain/seeds/postcing/postcing_3mm+tlrc.*  ./
