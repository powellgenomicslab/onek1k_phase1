##############################################################################
# Script information   

# Title: Generate conditioning plots
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "generate_conditioning_results_plot.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N cplots
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=12G
#$ -l tmp_requested=12G,tmpfree=12G
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

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# main script
RSCRIPT="/onek1k/independent_lead_eQTL_analysis/generate_conditining_results_plot.R"

# do the main job
$R_PATH/Rscript ${RSCRIPT} > rout/conditioning_plots.Rout


