##############################################################################
# Script information   

# Title: Generate locus plot data
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None

##############################################################################

# SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'
#$ -N smr
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

condition=$1
c

# Master sample file=
params="/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/${condition}_smr.lst"

celltype=`head -n $SGE_TASK_ID $params | tail -n 1 | awk '{print $1}'`
chr=`head -n $SGE_TASK_ID $params | tail -n 1 | awk '{print $2}'`
probe=`head -n $SGE_TASK_ID $params | tail -n 1 | awk '{print $3}'`

# SMR tool path
smr_tool="/onek1k/tools/smr_Linux"

# Plink data path
plink_dir="/onek1k/colocalisation_using_SMR/onek1kdata/plink_files_with_newIDs"

# GWAS data path
maindir="/onek1k/colocalisation_using_SMR/gwas_data"
gwas_path=${maindir}/${condition}.ma

# SMR result path
smr_path="/onek1k/colocalisation_using_SMR/smr_data"

${smr_tool} --bfile ${plink_dir}/plink_chr${chr}_renamed_201010 --gwas-summary ${gwas_path} \
--beqtl-summary ${smr_path}/${celltype}/${celltype}_chr${chr}_210211 --peqtl-smr 0.05 --diff-freq-prop 0.99  \
--out smr_analysis/smr_scatter_plot_data/${condition}_${celltype}_chr${chr} --plot --probe ${probe} --probe-wind 500 \
--gene-list smr_analysis/glist-hg19 