##############################################################################
# Script information   

# Title: Collate lead SNPs 
# Author: Seyhan Yazar
# Date: 2020-12-23
# Description: This bash script was written to collate lead SNPs from 
# conditional cis-eQTL mapping

##############################################################################

cd "/onek1k/single_cell_cis_eQTL_mapping"

mkdir All_Results

# lead SNPs
(head -1 Bmem/round1/Bmem_chr1_lead_eSNP.tsv && tail -n +2 -q Bmem/round1/Bmem_chr*_lead_eSNP.tsv ) > All_Results/Bmem_eSNP1.tsv
(head -1 CD4all/round1/CD4all_chr1_lead_eSNP.tsv && tail -n +2 -q CD4all/round1/CD4all_chr*_lead_eSNP.tsv ) > All_Results/CD4all_eSNP1.tsv
(head -1 CD4effCM/round1/CD4effCM_chr1_lead_eSNP.tsv && tail -n +2 -q CD4effCM/round1/CD4effCM_chr*_lead_eSNP.tsv ) > All_Results/CD4effCM_eSNP1.tsv
(head -1 CD4TGFbStim/round1/CD4TGFbStim_chr1_lead_eSNP.tsv && tail -n +2 -q CD4TGFbStim/round1/CD4TGFbStim_chr*_lead_eSNP.tsv ) > All_Results/CD4TGFbStim_eSNP1.tsv
(head -1 CD8all/round1/CD8all_chr1_lead_eSNP.tsv && tail -n +2 -q CD8all/round1/CD8all_chr*_lead_eSNP.tsv ) > All_Results/CD8all_eSNP1.tsv
(head -1 CD8eff/round1/CD8eff_chr1_lead_eSNP.tsv && tail -n +2 -q CD8eff/round1/CD8eff_chr*_lead_eSNP.tsv ) > All_Results/CD8eff_eSNP1.tsv
(head -1 CD8unknown/round1/CD8unknown_chr1_lead_eSNP.tsv && tail -n +2 -q CD8unknown/round1/CD8unknown_chr*_lead_eSNP.tsv ) > All_Results/CD8unknown_eSNP1.tsv
(head -1 DC/round1/DC_chr1_lead_eSNP.tsv && tail -n +2 -q DC/round1/DC_chr*_lead_eSNP.tsv ) > All_Results/DC_eSNP1.tsv
(head -1 MonoC/round1/MonoC_chr1_lead_eSNP.tsv && tail -n +2 -q MonoC/round1/MonoC_chr*_lead_eSNP.tsv ) > All_Results/MonoC_eSNP1.tsv
(head -1 MonoNC/round1/MonoNC_chr1_lead_eSNP.tsv && tail -n +2 -q MonoNC/round1/MonoNC_chr*_lead_eSNP.tsv ) > All_Results/MonoNC_eSNP1.tsv
(head -1 NKact/round1/NKact_chr1_lead_eSNP.tsv && tail -n +2 -q NKact/round1/NKact_chr*_lead_eSNP.tsv ) > All_Results/NKact_eSNP1.tsv
(head -1 NKmat/round1/NKmat_chr1_lead_eSNP.tsv && tail -n +2 -q NKmat/round1/NKmat_chr*_lead_eSNP.tsv ) > All_Results/NKmat_eSNP1.tsv
(head -1 Plasma/round1/Plasma_chr2_lead_eSNP.tsv && tail -n +2 -q Plasma/round1/Plasma_chr*_lead_eSNP.tsv ) > All_Results/Plasma_eSNP1.tsv
(head -1 BimmNaive/round1/BimmNaive_chr1_lead_eSNP.tsv && tail -n +2 -q BimmNaive/round1/BimmNaive_chr*_lead_eSNP.tsv ) > All_Results/BimmNaive_eSNP1.tsv
