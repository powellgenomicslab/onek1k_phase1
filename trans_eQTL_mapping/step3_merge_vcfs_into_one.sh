##############################################################################
# Script information   

# Title: Make a single vcf files for extract SNPs
# Author: Seyhan Yazar
# Date: 2021-03-02
# Description: None 

##############################################################################

## SGE SETTINGS
#$ -wd "/onek1k/trans_eQTL_mapping"
#$ -o "/onek1k/trans_eQTL_mapping/stdout"
#$ -S /bin/bash
#$ -q short.q
#$ -r yes
#$ -j no
#$ -l mem_requested=4G
#$ -l tmp_requested=4G
#$ -N merge_vcf

cd $SGE_0_WORKDIR
BIN=$SGE_O_WORKDIR/../../bin

conda activate /.conda/envs/onek1kEnv

cd /onek1k/trans_eQTL_mapping/chr_vcf_files

bcftools concat chr1_SNPs.recode.vcf \
    chr2_SNPs.recode.vcf \
    chr3_SNPs.recode.vcf \
    chr4_SNPs.recode.vcf \
    chr5_SNPs.recode.vcf \
    chr6_SNPs.recode.vcf \
    chr7_SNPs.recode.vcf \
    chr8_SNPs.recode.vcf \
    chr9_SNPs.recode.vcf \
    chr10_SNPs.recode.vcf \
    chr11_SNPs.recode.vcf \
    chr12_SNPs.recode.vcf \
    chr13_SNPs.recode.vcf \
    chr14_SNPs.recode.vcf \
    chr15_SNPs.recode.vcf \
    chr16_SNPs.recode.vcf \
    chr17_SNPs.recode.vcf \
    chr18_SNPs.recode.vcf \
    chr19_SNPs.recode.vcf \
    chr20_SNPs.recode.vcf \
    chr21_SNPs.recode.vcf \
    chr22_SNPs.recode.vcf -Oz > merged_cis_vcf_20210302.vcf.gz