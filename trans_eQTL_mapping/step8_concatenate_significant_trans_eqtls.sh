##############################################################################
# Script information   

# Title: Concatenate significant trans eQTLs
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "step8_concatenate_significant_trans_eqtls.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N sig_trans
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
RSCRIPT="/onek1k/trans_eQTL_mapping/step8_concatenate_significant_trans_eqtls.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} > rout/trans_combine.Rout
 