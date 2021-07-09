##############################################################################
# Script information   

# Title: Compare correlations between cell pairs
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for all cell pairs
# using "comparison_of_celltype_specific_correlations.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N comparison
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=24G
#$ -l tmp_requested=24G,tmpfree=24G
#$ -r yes
#$ -j y
#$ -wd '/onek1k/independent_lead_eQTL_analysis'
#$ -o '/onek1k/independent_lead_eQTL_analysis/stdout'

cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log
 
echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R upload
R_PATH="/.conda/envs/onek1kEnv/bin/"

# load cell type specific array.txt file
CTYPEFILE="/onek1k/independent_lead_eQTL_analysis/conditioning.batch2.txt" 

SAMPLE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
CELLTYPEA=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $2}'`
CELLTYPEB=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $3}'`

# log file
LOG=logs/log_validation_${CELLTYPE}

# main script
RSCRIPT="/onek1k/independent_lead_eQTL_analysis/comparison_of_celltype_specific_correlations.R"

# do the main job
$R_PATH/Rscript ${RSCRIPT} ${CELLTYPEA} ${CELLTYPEB} > rout/comparison_${CELLTYPEA}_${CELLTYPEB}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     echo -e `date` "\t${JOB_ID}\t${ARRAYID}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${ARRAYID}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS
