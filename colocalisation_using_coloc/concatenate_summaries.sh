##############################################################################
# Script information   

# Title: Concatenate summary associations 
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "concatenate_summaries.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N concat_sum
#$ -q short.q
#$ -l mem_requested=8G
#$ -l tmp_requested=8G,tmpfree=8G
#$ -l h='!epsilon*'
#$ -r yes
#$ -j yes
#$ -wd '/onek1k/colocalisation_using_coloc'
#$ -o '/onek1k/colocalisation_using_coloc/stdout'

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/disease.txt"

CONDITION=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# Create the output directory for the cell type
if [ ! -d coloc_analysis/summaries ]
then
    mkdir coloc_analysis/summaries
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/colocalisation_using_coloc/concatenate_summaries.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CONDITION} ${CELLTYPE} > rout/${CONDITION}_${CELLTYPE}.Rout
