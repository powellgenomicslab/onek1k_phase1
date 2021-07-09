##############################################################################
# Script information   

# Title: Run spearman's rank correlation tests for trans eqtl mapping
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: This bash script was written to run an array job for 14 cell types
# using "step6_run_spearman_corr_testing.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N trans
#$ -q short.q
#$ -l mem_requested=12G
#$ -l tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -wd '/onek1k/trans_eQTL_mapping'
#$ -o '/onek1k/trans_eQTL_mapping/stdout'

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Main script file
RSCRIPT="/onek1k/trans_eQTL_mapping/step6_run_spearman_corr_testing.R"

CELLTYPE=$1

if [ ! -d trans_eqtl_analysis/output2/${CELLTYPE} ]
then
    mkdir trans_eqtl_analysis/output2/${CELLTYPE}
else
     echo "Directory exists"
fi

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} > rout/trans_${CELLTYPE}.Rout
