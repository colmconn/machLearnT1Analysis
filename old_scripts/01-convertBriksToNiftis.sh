#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM
#set -x 

root=/data/sanDiego/

subjectsDataDir=$( readlink -f ../data )

subjects=$( cd $subjectsDataDir ; ls -1d *_A2 )

taskFile=convertBriksToNifti-TaskFile.txt
cat /dev/null > $taskFile
for subject in $subjects ; do

    echo "cd $subjectsDataDir/$subject; 3dcopy ${subject}+orig.HEAD ${subject}.anat.nii; 3dresample -orient RPI -inset ${subject}.anat.nii -prefix ${subject}.anat.std.nii" >> $taskFile
done


# jobname
# $ -N BrikToNii

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

    sge_command="qsub -N BrikToNii -q all.q -j y -m n -V -wd $( pwd ) -t 1-$nTasks" 
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
