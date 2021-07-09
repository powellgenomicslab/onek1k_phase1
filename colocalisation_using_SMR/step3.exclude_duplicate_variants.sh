##############################################################################
# Script information   

# Title: Remove duplicate variants
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "step3.exclude_duplicate_variants.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N smr_s3
#$ -q short.q
#$ -l mem_requested=12G
#$ -l tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -j y
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/$1"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# log file
LOG=logs/log_${CELLTYPE}

# Main script file
RSCRIPT="/onek1k/colocalisation_using_SMR/step3.exclude_duplicate_variants.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} ${CHRNUMBER} > rout/step3.exclude_dups_chr${CHRNUMBER}.Rout
