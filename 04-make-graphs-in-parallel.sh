#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM

## set -x 

root=/data/sanDiego/machLearnT1Analysis/data

subjects=$( cd $root; ls -1d [0-9][0-9][0-9]_A* )

CREATE_RSFC_GRAPHS=1
CREATE_ANAT_GRAPHS=0

[[ -d run ]] || mkdir run

taskFile="run/make-graphs-TaskFile.txt"
cat /dev/null > $taskFile
for subject in $subjects ; do

    ## 3dNetCorr partial correlation command line arg
    ## -part_corr
    ## 3dNetCorr fisher z command line arg    
    ## -fish_z
    
    if [[ $CREATE_RSFC_GRAPHS -eq 1 ]] ; then 
	preprocessed_rsfc_file=../data/$subject/rsfcPreprocessed/${subject}.pm.cleanEPI.MNI.nii.gz
	if [[ -f ${preprocessed_rsfc_file} ]] ; then 			\
	    echo "./make.netcorr.commands.r				\
	    	--source ${preprocessed_rsfc_file}			\
	    	--destination ../data/${subject}/rsfcGraphs		\
	    	--prefix ${subject}.pm.cleanEPI.aal2.whole.ts		\
	    	--window none						\
		--extra \\\\-part_corr					\
	    	--rois ../standard/aal2_for_SPM12/aal2.3mm.nii.gz	\
	    	--execute" >> $taskFile
	fi
    fi
	
    if [[ $CREATE_ANAT_GRAPHS -eq 1 ]] ; then 
	preprocessed_anat_file=../data/$subject/${subject}.anat_struc_GM_to_template_GM_mod_s2.nii.gz
	if [[ -f ${preprocessed_anat_file} ]] ; then 
	    echo "./make.netcorr.commands.r			\
	        --source ${preprocessed_anat_file}		\
        	--destination  ../data/${subject}/anatGraphs	\
	        --prefix  ${subject}.anat.aal2.whole		\
	        --window  none					\
	        --rois  ../standard/aal2_for_SPM12/aal2.nii.gz	\
        	--execute" >> $taskFile
	fi
    fi
done


# jobname
# $ -N makeGraphs

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

    sge_command="qsub -N makeGraphs -q all.q -j y -m n -V -wd $( pwd ) -o ../log/ -t 1-$nTasks" 
    #echo $sge_command
    echo "Queuing job... "
    ( exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF
)
    echo "Running qstat"
    qstat
    
fi
