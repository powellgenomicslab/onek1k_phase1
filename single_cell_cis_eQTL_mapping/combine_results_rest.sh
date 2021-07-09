##############################################################################
# Script information   

# Title: Collate eSNPs 
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This bash script was written to collate eSNP2 to 5 from 
# conditional cis-eQTL mapping

##############################################################################

cd "/onek1k/single_cell_cis_eQTL_mapping"

(head -1 Bmem/round$1/Bmem_chr2_eSNP$1.tsv && tail -n +2 -q Bmem/round$1/Bmem_chr*_eSNP$1.tsv ) > All_Results/Bmem_eSNP$1.tsv
(head -1 CD4all/round$1/CD4all_chr12_eSNP$1.tsv && tail -n +2 -q CD4all/round$1/CD4all_chr*_eSNP$1.tsv ) > All_Results/CD4all_eSNP$1.tsv
(head -1 CD4effCM/round$1/CD4effCM_chr2_eSNP$1.tsv && tail -n +2 -q CD4effCM/round$1/CD4effCM_chr*_eSNP$1.tsv ) > All_Results/CD4effCM_eSNP$1.tsv
(head -1 CD4TGFbStim/round$1/CD4TGFbStim_chr5_eSNP$1.tsv && tail -n +2 -q CD4TGFbStim/round$1/CD4TGFbStim_chr*_eSNP$1.tsv ) > All_Results/CD4TGFbStim_eSNP$1.tsv
(head -1 CD8all/round$1/CD8all_chr10_eSNP$1.tsv && tail -n +2 -q CD8all/round$1/CD8all_chr*_eSNP$1.tsv ) > All_Results/CD8all_eSNP$1.tsv
(head -1 CD8eff/round$1/CD8eff_chr3_eSNP$1.tsv && tail -n +2 -q CD8eff/round$1/CD8eff_chr*_eSNP$1.tsv ) > All_Results/CD8eff_eSNP$1.tsv
(head -1 CD8unknown/round$1/CD8unknown_chr2_eSNP$1.tsv && tail -n +2 -q CD8unknown/round$1/CD8unknown_chr*_eSNP$1.tsv ) > All_Results/CD8unknown_eSNP$1.tsv
(head -1 DC/round$1/DC_chr4_eSNP$1.tsv && tail -n +2 -q DC/round$1/DC_chr*_eSNP$1.tsv ) > All_Results/DC_eSNP$1.tsv
(head -1 MonoC/round$1/MonoC_chr6_eSNP$1.tsv && tail -n +2 -q MonoC/round$1/MonoC_chr*_eSNP$1.tsv ) > All_Results/MonoC_eSNP$1.tsv
(head -1 MonoNC/round$1/MonoNC_chr4_eSNP$1.tsv && tail -n +2 -q MonoNC/round$1/MonoNC_chr*_eSNP$1.tsv ) > All_Results/MonoNC_eSNP$1.tsv
(head -1 NKact/round$1/NKact_chr11_eSNP$1.tsv && tail -n +2 -q NKact/round$1/NKact_chr*_eSNP$1.tsv ) > All_Results/NKact_eSNP$1.tsv
(head -1 NKmat/round$1/NKmat_chr10_eSNP$1.tsv && tail -n +2 -q NKmat/round$1/NKmat_chr*_eSNP$1.tsv ) > All_Results/NKmat_eSNP$1.tsv
(head -1 Plasma/round$1/Plasma_chr3_eSNP$1.tsv && tail -n +2 -q Plasma/round$1/Plasma_chr*_eSNP$1.tsv ) > All_Results/Plasma_eSNP$1.tsv
(head -1 BimmNaive/round$1/BimmNaive_chr2_eSNP$1.tsv && tail -n +2 -q BimmNaive/round$1/BimmNaive_chr*_eSNP$1.tsv ) > All_Results/BimmNaive_eSNP$1.tsv