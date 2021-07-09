# Submit these qsub jobs to run association testing

cd "/onek1k/trans_eQTL_mapping"

qsub step6_run_spearman_corr_testing.sh BimmNaive
qsub step6_run_spearman_corr_testing.sh Bmem
qsub step6_run_spearman_corr_testing.sh CD4all
qsub step6_run_spearman_corr_testing.sh CD4effCM
qsub step6_run_spearman_corr_testing.sh CD4TGFbStim
qsub step6_run_spearman_corr_testing.sh CD8all
qsub step6_run_spearman_corr_testing.sh CD8eff
qsub step6_run_spearman_corr_testing.sh CD8unknown
qsub step6_run_spearman_corr_testing.sh DC
qsub step6_run_spearman_corr_testing.sh MonoC
qsub step6_run_spearman_corr_testing.sh MonoNC
qsub step6_run_spearman_corr_testing.sh NKact
qsub step6_run_spearman_corr_testing.sh NKmat
qsub step6_run_spearman_corr_testing.sh Plasma