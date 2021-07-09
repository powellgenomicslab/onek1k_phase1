##############################################################################
# Script information   

# Title: Conditional cis-eQTL mapping for CD4 NC cells - Round 5B
# Author: Seyhan Yazar
# Date: 2021-03-23
# Description: This bash script was written to run an array job for five subsets
# using "round5B.run_spearman_rank.sh" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/effect_of_cluster_size_on_eQTL_discovery'
#$ -o '/onek1k/effect_of_cluster_size_on_eQTL_discovery/stdout'
#$ -N power_round6
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
LOG=logs/log_round6

# set the variables
SAMPLE=$SGE_TASK_ID
PERCENT=$1

if [ ! -d effect_of_cluster_size_on_eQTL_discovery/percent${PERCENT}/round5B ]
then
     mkdir effect_of_cluster_size_on_eQTL_discovery/percent${PERCENT}/round5B
else
     echo "Directory exists"
fi
# Main script file
RSCRIPT="/onek1k/effect_of_cluster_size_on_eQTL_discovery/round5B.run_spearman_rank_test.R"

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} ${PERCENT} ${SAMPLE} > rout/round6_${PERCENT}_${SAMPLE}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS