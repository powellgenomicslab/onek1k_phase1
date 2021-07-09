##############################################################################
# Script information   

# Title: Identify significant trans eQTLs
# Author: Seyhan Yazar
# Date: 2021-03-05
# Description: This bash script was written to run an array job for 14 cell types
# using "step7_identify_significant_trans_eqtls.R" script. 

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
RSCRIPT="/onek1k/trans_eQTL_mapping/step7_identify_significant_trans_eqtls.R"

CELLTYPE=$1

if [ ! -d trans_eqtl_analysis/trans_eqtls_celltype_specific_output ]
then
    mkdir trans_eqtl_analysis/trans_eqtls_celltype_specific_output
else
     echo "Directory exists"
fi

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} > rout/trans_combine_${CELLTYPE}.Rout
