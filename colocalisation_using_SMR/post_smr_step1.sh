##############################################################################
# Script information   

# Title: Collate results and prepare plot data
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "post_smr_step1.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N pst_smr1
#$ -q short.q
#$ -l mem_requested=12G
#$ -l tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'

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
if [ ! -d smr_analysis/smr_results ]
then
    mkdir smr_analysis/smr_results
    mkdir smr_analysis/smr_results/smr_raw_results_condition_cell_level
    mkdir smr_analysis/smr_results/plot_summary_data_condition_cell_level
    mkdir smr_analysis/smr_results/plot_gene_list_condition_cell_level
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/colocalisation_using_SMR/post_smr_step1.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CONDITION} ${CELLTYPE} > rout/${CONDITION}_${CELLTYPE}.Rout
