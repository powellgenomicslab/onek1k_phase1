## SGE SETTINGS
#$ -S /bin/bash
#$ -N conditioning
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G
#$ -l tmp_requested=4G,tmpfree=4G
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

# R
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Master sample file
# Batch 1
# PARAMS="/onek1k/independent_lead_eQTL_analysis/conditioning.batch.txt"
# Batch 2
PARAMS="/onek1k/independent_lead_eQTL_analysis/conditioning.batch2.txt"

SAMPLE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE1=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`
CELLTYPE2=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $3}'`

CHRNUMBER=$1

# log file
LOG=logs/log_conditioning_${CELLTYPE1}_${CELLTYPE2}

# Create the output directory for the cell type
if [ ! -d expC1_snpC2_regression_results ]
then
   mkdir expC1_snpC2_regression_results
   mkdir residuals
   mkdir conditional_spearman_results
else
     echo "Directory exists"
fi


# main script
RSCRIPT="/onek1k/independent_lead_eQTL_analysis/conditioningBetweenCellTypes.R"

# do the main job
$R_PATH/Rscript ${RSCRIPT} ${CELLTYPE1} ${CELLTYPE2} ${CHRNUMBER} > rout/conditioning_${CELLTYPE1}_${CELLTYPE2}_$CHRNUMBER.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${SAMPLE}\t${CELLTYPE1}\t${CELLTYPE2}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${SAMPLE}\tconditioning\t${CELLTYPE}\t${CELLTYPE2}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS
