##############################################################################
# Script information   

# Title: Make .esi file for SMR analysis
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "make_bigesi_once.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'
#$ -N smr_s6
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Create the output directory for the cell type
if [ ! -d smr_analysis/smr_data/esi_files ]
then
     mkdir smr_analysis/smr_data/esi_files
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/colocalisation_using_SMR/make_bigesi_once.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${SGE_TASK_ID}
