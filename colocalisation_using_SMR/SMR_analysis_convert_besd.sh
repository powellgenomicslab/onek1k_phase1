##############################################################################
# Script information   

# Title: Convert .besd files to Matrix eQTL format
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "SMR_analysis_convert_besd.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/StdOut'
#$ -N smr
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=4G,tmp_requested=4G,tmpfree=4G
#$ -r yes
#$ -j y

CTYPEFILE="/onek1k/colocalisation_using_SMR/$1" 

SAMPLE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $2}'`

if [ ! -d onek1kdata/${CELLTYPE} ]
then
    mkdir onek1kdata/${CELLTYPE}
else
     echo "Directory exists"
fi

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# main script
RSCRIPT="/onek1k/colocalisation_using_SMR/SMR_analysis_convert_besd.R"

# do the main job
$R_PATH/Rscript ${RSCRIPT} ${CELLTYPE} ${SAMPLE} > Rout/SMR_analysis2_${CELLTYPE}_${SAMPLE}.Rout

# 
smr_tool="/onek1k/tools/smr_Linux"

cd onek1kdata/${CELLTYPE}

$smr_tool --eqtl-summary ${CELLTYPE}_chr${SAMPLE}_matrixeqtl_format.tsv --matrix-eqtl-format --make-besd --out ${CELLTYPE}_chr${SAMPLE}_besd
