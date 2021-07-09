##############################################################################
# Script information   

# Title: Prepare files for Matrix eQTL analysis
# Author: Seyhan Yazar
# Date: 2020-12-22
# Description: This bash script was written to run an array job for 14 cell types
# using "step1.prepare_files_for_MatrixEqtl.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -N smr_s1
#$ -q short.q
#$ -l mem_requested=32G
#$ -l tmp_requested=32G,tmpfree=32G
#$ -r yes
#$ -j y
#$ -wd '/onek1k/colocalisation_using_SMR'
#$ -o '/onek1k/colocalisation_using_SMR/stdout'

# debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

# Master sample file
PARAMS="/onek1k/Scripts/colocalisation_using_SMR/array_sets/$1"

CHRNUMBER=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $1}'`
CELLTYPE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{print $2}'`

# log file
LOG=logs/log_${CELLTYPE}

# Create the main directory for this analysis
if [ ! -d smr_analysis ]
then
     mkdir smr_analysis
     mkdir smr_analysis/matrix_eQTL
     mkdir smr_analysis/smr_data
else
     echo "Directory exists"
fi

# Create the output directory for the cell type
if [ ! -d smr_analysis/matrix_eQTL/${CELLTYPE} ]
then
     mkdir smr_analysis/matrix_eQTL/${CELLTYPE}
     mkdir smr_analysis/matrix_eQTL/${CELLTYPE}/input_files
else
     echo "Directory exists"
fi

# Main script file
RSCRIPT="/onek1k/colocalisation_using_SMR/step1.prepare_files_forMatrixEqtl.R"

# Do the main job
${R_PATH}/Rscript ${RSCRIPT} ${CELLTYPE} ${CHRNUMBER} > rout/step1.preparation_chr${CHRNUMBER}.Rout

# check exit status
STATUS=$?
if [[ $STATUS -eq 0 ]]; then
     # success, write MD5 verification file
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CHRNUMBER}\t${CELLTYPE}\tOK" >> $LOG
else
     echo -e `date` "\t${JOB_ID}\t${SGE_TASK_ID}\t${CHRNUMBER}\tarrayjob.one\t${CELLTYPE}\tFAIL\t${STATUS}" >> $LOG
fi
exit $STATUS