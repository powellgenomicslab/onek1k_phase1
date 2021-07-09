##############################################################################
# Script information   

# Title: Update epi files
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR/'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'
#$ -N smr_s6
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/$1"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# CHRNUMBER="22"
# CELLTYPE="BimmNaive"

smr_tool="/onek1k/tools/smr_Linux"
input_dir="/onek1k/colocalisation_using_SMR/smr_data/${CELLTYPE}"
epi_dir="/onek1k/colocalisation_using_SMR/smr_data/${CELLTYPE}/epi_files"

${smr_tool} --beqtl-summary ${input_dir}/${CELLTYPE}_chr${CHRNUMBER}_210211 --update-epi ${epi_dir}/${CELLTYPE}_chr${CHRNUMBER}.epi 