## SGE SETTINGS
#$ -S /bin/bash
#$ -wd '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/main_analysis'
#$ -o '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/main_analysis/stdout'
#$ -N supp_tbl
#$ -q short.q
#$ -pe smp 4
#$ -l mem_requested=12G,tmp_requested=12G,tmpfree=12G
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
R_PATH="/directflow/SCCGGroupShare/projects/SeyhanYazar/.conda/envs/onek1kEnv/bin/"

# Main script file
RSCRIPT="/directflow/SCCGGroupShare/projects/SeyhanYazar//onek1k/Scripts/abf3041_revision/main_analysis/all_eqtls_supplementary_table.R" 

# Do the main job
$R_PATH/Rscript --verbose ${RSCRIPT} > rout/supp_tbl.Rout
