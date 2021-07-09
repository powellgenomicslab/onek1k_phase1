##############################################################################
# Script information   

# Title: Make epi files
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "step5.make_bigepi.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout/'
#$ -N smr_s5
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/$1"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# Create the output directory for the cell type
if [ ! -d smr_analysis/smr_data/${CELLTYPE}/epi_files ]
then
     mkdir smr_analysis/smr_data/${CELLTYPE}/epi_files
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/colocalisation_using_SMR/step5.make_bigepi.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} ${CHRNUMBER}

