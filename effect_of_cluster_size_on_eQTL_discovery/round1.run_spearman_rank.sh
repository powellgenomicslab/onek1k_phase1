##############################################################################
# Script information   

# Title: Conditional cis-eQTL mapping for CD4 NC cells - Round 1
# Author: Seyhan Yazar
# Date: 2021-03-23
# Description: This bash script was written to run an array job for five subsets
# using "round1.run_spearman_rank.sh" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd 'effect_of_cluster_size_on_eQTL_discovery'
#$ -o 'effect_of_cluster_size_on_eQTL_discovery/stdout'
#$ -N power_round1
#$ -q short.q
#$ -pe smp 8
#$ -l mem_requested=12G,tmp_requested=12G,tmpfree=12G
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

# set the variables
SAMPLE=$SGE_TASK_ID
PERCENT=$1

if [ ! -d effect_of_cluster_size_on_eQTL_discovery/percent${PERCENT} ]
then
     mkdir effect_of_cluster_size_on_eQTL_discovery/percent${PERCENT}
     mkdir effect_of_cluster_size_on_eQTL_discovery/percent${PERCENT}/round1
     mkdir logs
     mkdir rout
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="effect_of_cluster_size_on_eQTL_discovery/onek1k/effect_of_cluster_size_on_eQTL_discovery/round1.run_spearman_rank.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${PERCENT} ${SAMPLE} > rout/round1_cd4${PERCENT}_${SAMPLE}.Rout

# log file
LOG=logs/log_round1

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${PERCENT}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS