#!/bin/bash

## set -x

studyName=machLearnT1Analysis

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/$studyName}
DATA=$ROOT/data
LOG_DIR=$ROOT/log
SCRIPTS_DIR=${ROOT}/scripts

. ${SCRIPTS_DIR}/logger_functions.sh

if [[ $# -gt 0 ]] ; then
    subjects="$*"
else

    subjects="130_A 136_A 140_A 149_A 161_A 165_A
167_A 300_A 301_A 307_A 309_A 311_A 316_A 324_A 326_A 328_A 329_A
336_A 339_A 364_A 368_A"
    
    ## subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.nat.txt )"
    subjects=$( cd ../data; ls -1d *_A2 )
    # subjectCount=$( cd ../data; ls -1d *_A2 | wc -l )
fi

taskFile=$SCRIPTS_DIR/run/convertDicoms-TaskFile.$BASHPID

cat /dev/null > $taskFile


jobname="convertDicoms"
(( i=1 ))
for subject in ${subjects} ; do
    info_message "$( printf "Adding script(s) for subject %s (%03d of %03d) to task file\n" $subject $i $subjectCount )"

    echo "$SCRIPTS_DIR/00-convertDicoms.sh -s $subject" >> ${taskFile}

    ## echo "./03-singleSubjectRsfc.sh -s $subject -l ../data/config/dlpfc.seed.list.txt" >> ${taskFile}
    ## echo "./03-singleSubjectRsfc.sh -s $subject -l ../data/config/seed.list.txt" >> ${taskFile}
    (( i=i+1 ))
done

## jobname
#$ -N $jobname

## queue
#$ -q all.q

## binary? 
#$ -b y

## rerunnable?
#$ -r y

## merge stdout and stderr?
#$ -j y

## send no mail
#$ -m n

## execute from the current working directory
#$ -cwd

## use a shell to run the command
#$ -shell yes 

## set the shell
#$ -S /bin/bash

## preserve environment
#$ -V 

[[ ! -d $LOG_DIR ]] && mkdir $LOG_DIR

nTasks=$( cat $taskFile | wc -l )
sge_command="qsub -N $jobname -q all.q -j y -m n -V -wd $( pwd ) -o $LOG_DIR -t 1-$nTasks" 
echo $sge_command
( exec $sge_command <<EOF
#!/bin/sh

#$ -S /bin/sh

command=\`sed -n -e "\${SGE_TASK_ID}p" $taskFile\`

exec /bin/sh -c "\$command"
EOF
)

echo "Running qstat"
qstat
