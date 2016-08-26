#!/bin/bash

set -x 

#$ -S /bin/bash

export PYTHONPATH=/data/software/afni

## use the newer faster despiking method. comment this out to get the
## old one back
export AFNI_3dDespike_NEW=YES

## only use a single thread since we're going to run so many subjects
## in parallel
export OMP_NUM_THREADS=8

subject="$@"

cd /data/sanDiego/machLearnT1Analysis/data/${subject}

outputDir=alignmentTest
rm -rf alignmentTest

mkdir alignmentTest
cd alignmentTest/
cp ../${subject}.resting+orig.* ./
if [[ -f ../${subject}.anat_clp+orig.HEAD ]] ; then
    3dcopy ../${subject}.anat_clp+orig. ${subject}.anat_clp
    anatFile=${subject}.anat_clp+orig.HEAD
else
    cp ../${subject}.anat+orig.* ./
    anatFile=${subject}.anat+orig.HEAD
fi

## cp ../${subject}.anat+orig.* ./
## orientation=$( 3dinfo -orient ${subject}.anat.nii.gz )
# if [[ "${orientation}" != "RPI" ]] ; then
#     echo "*** Reorienting anatomy to RPI"
#     3dresample -orient RPI -inset ${subject}.anat.nii.gz -prefix ${subject}.anat.std.nii
#     anatFile=${subject}.anat.std.nii.gz
# else
#     anatFile=${subject}.anat.nii.gz
# fi

## anatFile=${subject}.anat+orig.HEAD
3dZeropad -I 30 -S 30 -prefix ${anatFile%%+*}.zp ${anatFile}
anatFile=${anatFile%%+*}.zp+orig.HEAD						      

## 3dTcat -prefix ${subject}.resting.tcat ${subject}.resting+orig.'[3..$]'
3dcopy ${subject}.resting+orig ${subject}.resting.tcat
3dDespike -NEW -nomask -prefix ${subject}.resting.despike ${subject}.resting.tcat+orig.
3dTshift -tzero 0 -quintic -prefix ${subject}.resting.tshift ${subject}.resting.despike+orig.

3dZeropad -I 30 -S 30 -prefix ${subject}.resting.tshift.zp ${subject}.resting.tshift+orig
epiFile=${subject}.resting.tshift.zp+orig.HEAD						      

align_epi_anat.py -anat2epi			\
		  -anat ${anatFile}		\
		  -epi ${epiFile}		\
		  -epi_base 0			\
		  -volreg off			\
		  -tshift off			\
		  -cost lpc			\
		  -giant_move			\
		  -multi_cost lpa lpc+ZZ mi





# afni_proc.py -subj_id ${subject}								\
#              -script ${preprocessingScript}							\
# 	     -out_dir ${outputDir}								\
# 	     -scr_overwrite									\
# 	     -blocks despike tshift align tlrc volreg						\
# 	     -copy_anat /data/sanDiego/rsfcGraphAnalysis/data/${subject}/${subject}.anat.nii.gz	\
# 	     -dsets /data/sanDiego/rsfcGraphAnalysis/data/${subject}/${subject}.resting+orig.HEAD	\
# 	     -tcat_remove_first_trs 3								\
# 	     -tlrc_base MNI_caez_N27+tlrc							\
# 	     -volreg_align_to MIN_OUTLIER							\
# 	     -volreg_tlrc_warp									\
# 	     -execute

