##############################################################################
# Script information   

# Title: Combine conditioning results
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for cell pairs
# using "combine_conditioning_results.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/independent_lead_eQTL_analysis'
#$ -o '/onek1k/independent_lead_eQTL_analysis/stdout'
#$ -N combine_conditioning
#$ -q short.q
#$ -pe smp 2
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
LOG=logs/log_combining_results

# load cell type specific array.txt file
CONDFILE="/onek1k/independent_lead_eQTL_analysis/conditioning.batch.txt" 

SAMPLE=`head -n $SGE_TASK_ID $CONDFILE | tail -n 1 | awk '{print $1}'`
CELLLABEL1=`head -n $SGE_TASK_ID $CONDFILE | tail -n 1 | awk '{print $2}'`
CELLLABEL2=`head -n $SGE_TASK_ID $CONDFILE | tail -n 1 | awk '{print $3}'`

if [ ! -d cs_all_results ]
then
     mkdir cs_all_results
else
     echo "Directory exists"
fi


# main script
RSCRIPT="/onek1k/independent_lead_eQTL_analysis/combine_conditioning_results.R"

# do the main job
$R_PATH/Rscript ${RSCRIPT} ${CELLLABEL1} ${CELLLABEL2} > rout/combine_${CELLLABEL1}_${CELLLABEL2}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL1}\t${CELLLABEL2}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CELLLABEL1}\t${CELLLABEL2}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS
