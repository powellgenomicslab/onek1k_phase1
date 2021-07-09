##############################################################################
# Script information   

# Title: Convert vcf file to plink format
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: None

##############################################################################

## SGE SETTINGS
#$ -wd "/onek1k/trans_eqtl_mapping"
#$ -o "/onek1k/trans_eqtl_mapping/stdout"
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -l mem_requested=4G
#$ -l tmp_requested=4G
#$ -N convert

cd $SGE_0_WORKDIR
BIN=$SGE_O_WORKDIR/../../bin

conda activate /.conda/envs/onek1kEnv

vcfFile="/onek1k/trans_eqtl_mapping/chr_vcf_files/merged_cis_vcf_20210302.vcf.gz"
outDir="onek1k/trans_eqtl_mapping

plink2 --vcf ${vcfFile} --const-fid --make-bed --out ${outDir}/snps_20210302