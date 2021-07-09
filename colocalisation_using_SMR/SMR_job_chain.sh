##############################################################################
# Script information   

# Title: Prepare SMR files
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: Seven steps of SMR analysis preparations are chained to run as a single job 

##############################################################################

## SGE SETTINGS
#$ -S /bin/bash
#$ -cwd
#$ -N jobChain
#$ -r yes
#$ -j y

scripts_directory="/onek1k/colocalisation_using_SMR"

one=$(qsub -t 19:19 $scripts_directory/step1.submit.arrayjob.sh $1)
echo $one
two=$(qsub -t 19:19 -hold_jid smr_s1 $scripts_directory/step2.submit.arrayjob.sh $1)
echo $two
three=$(qsub -t 19:19 -hold_jid smr_s2 $scripts_directory/step3.exclude_duplicate_variants.sh $1)
echo $three
four=$(qsub -t 19:19 -hold_jid smr_s3 $scripts_directory/step4.make_besd.sh $1)
echo $four
five=$(qsub -t 19:19 -hold_jid smr_s4 $scripts_directory/step5.make_bigepi.sh $1)
echo $five
six=$(qsub -t 19:19 -hold_jid smr_s5 $scripts_directory/step6.update_epi.sh $1)
echo $six
seven=$(qsub -t 19:19 -hold_jid smr_s6 $scripts_directory/step7.update_esi.sh $1)
echo $seven