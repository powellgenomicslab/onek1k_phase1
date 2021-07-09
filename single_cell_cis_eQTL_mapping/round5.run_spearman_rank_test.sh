##############################################################################
# Script information   

# Title: Conditional cis-eQTL mapping - round 5
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This bash script was written to run an array job 
# per chromosome for 14 cell types using "round5.run_spearman_rank_test.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/single_cell_cis_eQTL_mapping'
#$ -o '/onek1k/single_cell_cis_eQTL_mapping/stdout'
#$ -N j21_round5
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log
 
echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# log file
LOG=logs/log_round5

# load cell type specific array.txt file
CTYPEFILE="/onek1k/Scripts/abf3041_revision/main_analysis/$1"

SAMPLE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $2}'`

if [ ! -d ${CELLLABEL}/round5 ]
then
     mkdir ${CELLLABEL}/round5
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/Scripts/single_cell_cis_eQTL_mapping/round5.run_spearman_rank_test.R" 

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} ${CELLLABEL} ${SAMPLE} > rout/round5_${CELLLABEL}_${SAMPLE}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS