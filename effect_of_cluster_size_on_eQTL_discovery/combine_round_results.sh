cd "/onek1k/effect_of_cluster_size_on_eQTL_discovery"

# mkdir percent_specific_results

# (head -1 percent1/round1/percent1_chr1_lead_eSNP.tsv && tail -n +2 -q percent1/round1/percent1_chr*_lead_eSNP.tsv ) > percent_specific_results/percent1_eSNP1.tsv
# (head -1 percent5/round1/percent5_chr1_lead_eSNP.tsv && tail -n +2 -q percent5/round1/percent5_chr*_lead_eSNP.tsv ) > percent_specific_results/percent5_eSNP1.tsv
# (head -1 percent10/round1/percent10_chr1_lead_eSNP.tsv && tail -n +2 -q percent10/round1/percent10_chr*_lead_eSNP.tsv ) > percent_specific_results/percent10_eSNP1.tsv
# (head -1 percent25/round1/percent25_chr1_lead_eSNP.tsv && tail -n +2 -q percent25/round1/percent25_chr*_lead_eSNP.tsv ) > percent_specific_results/percent25_eSNP1.tsv
# (head -1 percent50/round1/percent50_chr1_lead_eSNP.tsv && tail -n +2 -q percent50/round1/percent50_chr*_lead_eSNP.tsv ) > percent_specific_results/percent50_eSNP1.tsv
# (head -1 percent75/round1/percent75_chr1_lead_eSNP.tsv && tail -n +2 -q percent75/round1/percent75_chr*_lead_eSNP.tsv ) > percent_specific_results/percent75_eSNP1.tsv

(head -1 percent1/round$1/percent1_chr2_eSNP$1.tsv && tail -n +2 -q percent1/round$1/percent1_chr*_eSNP$1.tsv ) > percent_specific_results/percent1_eSNP$1.tsv
(head -1 percent5/round$1/percent5_chr2_eSNP$1.tsv && tail -n +2 -q percent5/round$1/percent5_chr*_eSNP$1.tsv ) > percent_specific_results/percent5_eSNP$1.tsv
(head -1 percent10/round$1/percent10_chr2_eSNP$1.tsv && tail -n +2 -q percent10/round$1/percent10_chr*_eSNP$1.tsv ) > percent_specific_results/percent10_eSNP$1.tsv
(head -1 percent25/round$1/percent25_chr2_eSNP$1.tsv && tail -n +2 -q percent25/round$1/percent25_chr*_eSNP$1.tsv ) > percent_specific_results/percent25_eSNP$1.tsv
(head -1 percent50/round$1/percent50_chr2_eSNP$1.tsv && tail -n +2 -q percent50/round$1/percent50_chr*_eSNP$1.tsv ) > percent_specific_results/percent50_eSNP$1.tsv
(head -1 percent75/round$1/percent75_chr2_eSNP$1.tsv && tail -n +2 -q percent75/round$1/percent75_chr*_eSNP$1.tsv ) > percent_specific_results/percent75_eSNP$1.tsv