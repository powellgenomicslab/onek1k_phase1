##############################################################################
# Script information   

# Title: Make besd file for SMR script
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'
#$ -N smr_s4
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=12G,tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -j y

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/$1"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# CHRNUMBER="10"
# CELLTYPE="CD4all"

smr_tool="/onek1k/tools/smr_Linux"
input_dir="/onek1k/colocalisation_using_SMR/matrix_eQTL/${CELLTYPE}/output_files"
output_dir="colocalisation_using_SMR/smr_data/${CELLTYPE}"

mkdir $output_dir

$smr_tool --eqtl-summary ${input_dir}/${CELLTYPE}_chr${CHRNUMBER}_cis_eqtls_210211_uniq.tsv --matrix-eqtl-format --make-besd --out ${output_dir}/${CELLTYPE}_chr${CHRNUMBER}_210211