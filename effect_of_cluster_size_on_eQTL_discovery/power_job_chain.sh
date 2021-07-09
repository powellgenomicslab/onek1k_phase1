# Submit qsub jobs in a chain 

## SGE SETTINGS
#$ -S /bin/bash
#$ -cwd
#$ -N jobChain
#$ -r yes
#$ -j y

scripts_directory="/onek1k/effect_of_cluster_size_on_eQTL_discovery"

one=$(qsub -t 6:6 $scripts_directory/round2.run_spearman_rank_test.sh $1)
echo $one
two=$(qsub -t 6:6 -hold_jid power_round2 $scripts_directory/round3.run_spearman_rank_test.sh $1)
echo $two
three=$(qsub -t 6:6 -hold_jid power_round3 $scripts_directory/round4.run_spearman_rank_test.sh $1)
echo $three
four=$(qsub -t 6:6 -hold_jid power_round4 $scripts_directory/round5.run_spearman_rank_test.sh $1)
echo $four
five=$(qsub -t 6:6 -hold_jid power_round5 $scripts_directory/round5B.run_spearman_rank_test.sh $1)
echo $five
