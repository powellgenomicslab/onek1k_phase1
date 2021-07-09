# Run the following scripts to generate locus plot data

### AS ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/as_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh as

### CROHNS ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/crohns_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh crohns

### IBD ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/ibd_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh ibd

### MS ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/ms_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh ms

### RA ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/ra_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh ra

### SLE ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/sle_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh sle

### T1DM ###
file='/onek1k/colocalisation_using_SMR/smr_results/plot_gene_list_condition_level/t1dm_smr.lst'
num_of_lines=$(cat $file | wc -l)
qsub -t 1:$num_of_lines post_smr_step3_locus_plot_data_generation.sh t1dm