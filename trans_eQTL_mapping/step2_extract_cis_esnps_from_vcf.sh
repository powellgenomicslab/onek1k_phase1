##############################################################################
# Script information   

# Title: Extract cis SNPs from the vcf files
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
#$ -l mem_requested=40G
#$ -l tmp_requested=40G
#$ -N snp.extraction

cd $SGE_0_WORKDIR
BIN=$SGE_O_WORKDIR/../../bin
export PATH=$PATH:$BIN
export PERL5LIB=/share/ClusterShare/software/contrib/sccg/vcftools/src/perl
export PATH=$PATH:/share/ClusterShare/software/contrib/sccg/vcftools/bin

vcfDir="/onek1k/imputed_data/filter_vcf_r08_maf005"
outputDir="/onek1k/science_revision_Dec20/trans_eqtl_mapping/chr_vcf_files"
snpFile="/onek1k/science_revision_Dec20/trans_eqtl_mapping/OneK1K_cis_eSNP_list.txt"

# Create the output directory for the cell type
if [ ! -d ${outputDir} ]
then
     mkdir ${outputDir}
else
     echo "Directory exists"
fi

# Extract eSNPs from each file

for CHR in `seq 1 22`;
do
vcftools --gzvcf ${vcfDir}/chr${CHR}.dose.filtered.R2_08_maf005.vcf.gz --snps ${snpFile} --recode --recode-INFO-all \
    --out ${outputDir}/chr${CHR}_SNPs;
done