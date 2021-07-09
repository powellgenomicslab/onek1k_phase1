##############################################################################
# Script information   

# Title: Identify PEER factors
# Author: Seyhan Yazar
# Date: 2020-12-22
# Description: This bash script was written to run an array job for 14 cell types
# using "identify_peer_factors.R" script. 
# Note: PEER package is only available for earlier versions of R. Please generate
# a separate environment for it. 

##############################################################################

## SGE SETTINGS ##
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/stdout
#$ -N peer_test
#$ -q short.q
#$ -pe smp 8
#$ -l mem_requested=12G,tmp_requested=12G,tmpfree=12G
#$ -r yes
#$ -j y

# Debug
set -x

# Clear the environment
. /etc/profile.d/modules.sh

echo $PATH

# Activate conda environment where PEER package is installed
conda activate peerEnv

# Set the text file for array job
CTYPEFILE='celltypes.txt'

SAMPLE=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $1}'`
CELLLABEL=`head -n $SGE_TASK_ID $CTYPEFILE | tail -n 1 | awk '{print $2}'`

# Main script file
RSCRIPT="identify_PEER_factors.R"

# Run the main script
Rscript ${RSCRIPT} ${CELLLABEL}