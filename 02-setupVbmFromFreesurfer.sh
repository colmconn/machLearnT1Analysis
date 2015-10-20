#!/bin/bash

# set -x

SUBJECTS_DIR=/data/sanDiego/freesurferAnalysis/data

studyRoot=/data/sanDiego/machLearnT1Analysis
dataRoot=$studyRoot/data/
vbmDir=$dataRoot/vbm

if [[ ! -d $vbmDir/struc ]] ;then
    mkdir -p $vbmDir/struc
fi

if [[ ! -d ../log ]] ; then 
    mkdir ../log 
fi

subjects=$( cd ../data ; ls -d1 *_A2 )
#subjects="105_A"
echo $subjects

function makeStdAnatLink() {

    local subject="$1"

    stdAnatFile=$dataRoot/${subject}/${subject}.anat.std.nii.gz
    if [[ ! -f $stdAnatFile ]] ; then
	echo "Couldn't make link for $subject: No such file $stdAnatFile"
    else
	( cd $vbmDir;  ln -sf $stdAnatFile        ${subject}.anat.nii.gz ; symlinks -c ./ )
    fi
}

taskFile=convertMgzToNii-TaskFile.txt
cat /dev/null > $taskFile
for subject in ${subjects} ; do
    
    if [[ ! -f $SUBJECTS_DIR/${subject}/mri/brainmask.mgz ]] ; then
	problemBrains="$subject $problemBrains"
    else
	makeStdAnatLink $subject
	
	echo "( cd $vbmDir/struc; \
	  mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/T1.mgz ${subject}.anat_struc.RSP.nii; \
	  3dresample -orient RPI -inset ${subject}.anat_struc.RSP.nii -prefix ${subject}.anat_struc.nii ; \
	  mri_convert -it mgz -ot nii $SUBJECTS_DIR/${subject}/mri/brainmask.mgz  ${subject}.anat_struc_brain.RSP.nii ; \
	  3dresample -orient RPI -inset ${subject}.anat_struc_brain.RSP.nii -prefix ${subject}.anat_struc_brain.nii ; \
	  rm -f ${subject}.anat_struc.RSP.nii ; \
	  rm -f ${subject}.anat_struc_brain.RSP.nii )" >> $taskFile
    fi
done
		    
echo "Print problem brains: $problemBrains"

# jobname
# $ -N MgzToNii

# queue
# $ -q all.q

# binary?
# $ -b y

# rerunnable?
# $ -r y

# merge stdout and stderr?
# $ -j y

# send no mail
# $ -m n

# execute from the current working directory
# $ -cwd

# use a shell to run the command
# $ -shell yes 
# set the shell
# $ -S /bin/bash

# preserve environment
# $ -V 

if [[ 1 == 1 ]] ; then
    nTasks=$( cat $taskFile | wc -l )

    sge_command="qsub -o ../log/$MgzToNii.log -N MgzToNii -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
    #echo $sge_command
    echo "Queuing job... "
    exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF

    echo "Running qstat"
    qstat

fi
