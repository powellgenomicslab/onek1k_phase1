##############################################################################
# Script information   

# Title: Run cell composition GWAS
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd 'onek1k/single_cell_cis_eQTL_mapping '
#$ -o 'onek1k/single_cell_cis_eQTL_mapping stdout'
#$ -N cell_prop_gwas_plink
#$ -q short.q
#$ -l mem_requested=32G,tmp_requested=32G,tmpfree=32G
#$ -r yes
#$ -j y

cat /tmp/prolog_exec_"$JOB_ID"_"$SGE_TASK_ID".log
 
echo "JOB: $JOB_ID TASK: $SGE_TASK_ID"
echo "$HOSTNAME $tmp_requested $TMPDIR"

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# Upload modules
module load /share/ClusterShare/Modules/modulefiles/contrib/evaben/gcc/gcc-7.3.0/7.3.0
PLINK_PATH="/onek1k/tools"

# log file
LOG=Log_Files/log_cell_prop_plink_gwas

# load cell type specific array.txt file
CTYPEFILE="/onek1k/single_cell_cis_eQTL_mapping/$1" 

SAMPLE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $2}'`

# Phenotype file
PHENO="/onek1k/Plink_Phenotypes/${CELLLABEL}_phenotype.tsv"

# Covariate file
COVAR="/onek1k/single_cell_cis_eQTL_mapping/Plink_Covariates/${CELLLABEL}_covariates.txt"

# Genotype file
GENO="/onek1k/imputed_data/filter_vcf_r08_maf005/plink_chr${SAMPLE}"

# Results Folder
if [ ! -d /onek1k/single_cell_cis_eQTL_mapping/Plink_Outputs/${CELLLABEL} ]
then
     mkdir /onek1k/single_cell_cis_eQTL_mapping/Plink_Outputs/${CELLLABEL}
else
     echo "Directory exists"
fi

# Run plink
$PLINK_PATH/plink --bfile ${GENO} --pheno ${PHENO} --covar ${COVAR} --allow-no-sex \
    --linear hide-covar --out /onek1k/single_cell_cis_eQTL_mapping/Plink_Outputs/${CELLLABEL}/${CELLLABEL}_$SAMPLE