##############################################################################
# Script information   

# Title: Calculate average expression per person
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: This bash script was written to run an array job for 14 cell types
# using "step3.exclude_duplicate_variants.R" script. 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/onek1k/effect_of_cluster_size_on_eQTL_discovery'
#$ -o '/onek1k/effect_of_cluster_size_on_eQTL_discovery/stdout'
#$ -N calc_means
#$ -q short.q
#$ -l mem_requested=196G
#$ -l tmp_requested=196G,tmpfree=196G
#$ -r yes
#$ -j y

# debug
set -x

# Clear environment
. /etc/profile.d/modules.sh

# R path
R_PATH="/.conda/envs/onek1kEnv/bin/"

PERCENT=$1

if [ ! -d effect_of_cluster_size_on_eQTL_discovery/count_averages ]
then
     mkdir effect_of_cluster_size_on_eQTL_discovery/count_averages
else
     echo "Directory exists"
fi

# Main script
RSCRIPT="/onek1k/effect_of_cluster_size_on_eQTL_discovery/calculate_average_expression_per_person.R"

${R_PATH}/Rscript ${RSCRIPT} ${PERCENT} > rout/average.rout