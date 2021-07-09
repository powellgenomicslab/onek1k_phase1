##############################################################################
# Script information   

# Title: Run SMR analysis
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR'
#$ -N smr
#$ -q short.q
#$ -pe smp 16
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y
# -l h='!test*'

# Master sample file
PARAMS="/onek1k/colocalisation_using_SMR/array_sets/disease.txt"

CONDITION=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

smr_tool="/onek1k/tools/smr_Linux"
plink_dir="/onek1k/colocalisation_using_SMR/onek1kdata/plink_files_with_newIDs"
gwas_file_path="/onek1k/colocalisation_using_SMR/gwas_data/${CONDITION}.ma"
smr_path="/onek1k/colocalisation_using_SMR/smr_data"

# Create the output directory for the cell type
if [ ! -d smr_analysis/smr_output/${CONDITION} ]
then
      mkdir smr_analysis/smr_output/${CONDITION}
else
     echo "Directory exists"
fi

mkdir smr_analysis/smr_output/${CONDITION}/${CELLTYPE}

cd /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/smr_output/${CONDITION}/${CELLTYPE}

for chr in {1..22}; do
     ${smr_tool} --bfile ${plink_dir}/plink_chr${chr}_renamed_201010 --gwas-summary ${gwas_file_path} \
      --beqtl-summary ${smr_path}/${CELLTYPE}/${CELLTYPE}_chr${chr}_210211 --peqtl-smr 0.05 \
      --out ${CELLTYPE}_chr${chr} --thread-num 16 --diff-freq-prop 0.99
done
