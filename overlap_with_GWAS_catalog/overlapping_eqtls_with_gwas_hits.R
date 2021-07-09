##############################################################################
# Script information                                                      
# Title: Overlapping eQTLs with GWAS hits
# Author: Seyhan Yazar
# Date: 2021-03-12
# Description: Includes supplementary figure S22
##############################################################################

# MUST HAVES
# Import libraries
Packages <- c("plyr", "dplyr", "ggplot2", "reshape2", "magrittr",
              "purrr", "tidyr", "UpSetR", "RColorBrewer", "cowplot",
              "gridExtra", "ggpubr", "ComplexHeatmap", "tibble",
              "tidyverse", "vctrs", "data.table", "ggthemes", "likert", 
              "epiR", "zoo", "ggpmisc")
invisible(lapply(Packages, library, character.only = TRUE))

# Colour palettee
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
# Order of independent eQTLs
rank_order <- c("eSNP1","eSNP2","eSNP3", "eSNP4","eSNP5")

# Order of cells
cell_order <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")

cell_labels <- c("CD4 Naïve&CM","CD4 EM&TEMRA","CD4 SOX4","CD8 EM&TEMRA","CD8 Naïve&CM","CD8 S100B","NK", "NK Recruiting","Plasma",
                 "B Mem","B Imm&Naïve","Mono C","Mono NC", "DC")

# Set order of plotting
corder <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")

# snp hrc ids
df_hrc <- fread("/onek1k/hrc_ids_eqtl.tsv")
colnames(df_hrc)

# Folder to save figures nd tables
FigTabDir <- "/onek1k/Figures_and_Tables"
# dir.create(FigTabDir)

cells <- data.frame(cell_order,cell_labels,corder)
colnames(cells)[1] <- "celltype"
cells$celltype <- as.character(cells$celltype)

setwd("/onek1k/All_Results")
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { tryCatch(read_tsv(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames

# Add hrc ids to all list
dataset <- lapply(dataset, function(x) left_join(x,df_hrc,by="snpid"))

# Combine lists
dataset2<- bind_rows(dataset, .id="list_name")
dataset2$celltype <- gsub("_.*","", dataset2$list_name)
dataset2$eSNP_rank <- gsub(".*_","",dataset2$list_name)
dataset2$celltype <- factor(dataset2$celltype, levels=cell_order)
dataset2 <- left_join(dataset2, cells, by="celltype")

df_eqtl <- dataset2 %>% select(corder, geneid, ID, snpid, estimate, p.value, localFDR, eSNP_rank)
colnames(df_eqtl) <- c("cell type")

# Add these rsIDs 
dataset2[dataset2$snpid=="17:44355634_A", "ID"] <- "rs2696531"

gwas <- fread("/onek1k/GWAS_Results/gwas_results_with_hrcids_20200212.tsv")

# Function to identify overlapping genes by disease/trait column
overlap <-  function (x, gwas, analysis) {
  disease <- gwas[(gwas$'DISEASE/TRAIT'==x),]
  
  disease <- disease[,c(12,13,14,15,21,37)]
  disease_genes <- unique(unlist(strsplit(disease$`REPORTED GENE(S)`, ", ")))
  no <- length(disease_genes)
  study <- unique(disease$`STUDY ACCESSION`)
  nostudy <- length(study)
  disease_df <- analysis[(analysis$geneid %in% disease_genes),]
  nogenes <- disease_df %>% group_by(corder) %>% count()
  nogenes$tot <- no
  nogenes$perc <- (nogenes$n/nogenes$tot)*100
  nogenes$study <- nostudy 
  nogenes
}

# Apply overlap function to 8 immune diseases
ms <- overlap("Multiple sclerosis", gwas, dataset2) %>% data.frame()
diabetes1 <- overlap("Type 1 diabetes", gwas, dataset2) %>% data.frame()
# diabetes2 <- overlap("Type 2 diabetes", gwas,dataset2) %>% data.frame()
crohn <- overlap("Crohn's disease", gwas, dataset2) %>% data.frame()
sle <- overlap("Systemic lupus erythematosus",gwas, dataset2) %>% data.frame()
ibd <- overlap("Inflammatory bowel disease", gwas, dataset2) %>% data.frame()
ra <- overlap("Rheumatoid arthritis", gwas, dataset2) %>% data.frame()
as <- overlap("Ankylosing spondylitis", gwas, dataset2) %>% data.frame()

diseaseList <- list(MS=ms,T1DM=diabetes1,CD=crohn,SLE=sle,IBD=ibd,RA=ra,AS=as)

order_to_implement <- c("CD4 all","CD4 Eff&CM", "CD4 TGFbstim" , "CD8 Eff", "CD8 all", "CD8 other",
                        "NK mat", "NK act","Plasma", "B Mem", "B imm&Naive", "Mono c", "Mono nc", "DC")

diseaseList <- lapply (diseaseList, function(x) x %>% arrange(factor(corder)))

disease_merged <- bind_rows(diseaseList, .id="disease")
disease_merged$level_order <- factor(disease_merged$corder, levels=corder)
disease_merged <- disease_merged %>% arrange(corder)

# Set a number of 'empty bar'
empty_bar <- 3

# Add lines to the initial dataset
to_add <- data.frame(matrix(NA, empty_bar*nlevels(factor(disease_merged$disease)), ncol(disease_merged)))
colnames(to_add) <- colnames(disease_merged)
to_add$disease <- rep(levels(factor(disease_merged$disease)), each=empty_bar)
disease_merged <- rbind(disease_merged, to_add)
disease_merged <- disease_merged %>% arrange(disease)
disease_merged$id <- seq(1, nrow(disease_merged))

# Get the name and the y position of each label
label_data <- disease_merged
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- disease_merged %>% 
  group_by(disease) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
# grid_data <- grid_data[-1,]

# Make the plot
p_genes <- ggplot(disease_merged, aes(x=as.factor(id), y=perc, color=disease)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill="white", color="black") +
  # geom_bar(aes(x=as.factor(id), y=perc), stat="identity", alpha=0.5) +
  # scale_fill_manual(values=c("#000000", "#000000","#000000", "#000000", "#000000", "#000000","#000000")) +
  #scale_color_manual(values=c("#999999", "#999999", "#999999","#999999","#999999","#999999","#999999")) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(disease_merged$id),4), y = c(10,20, 30, 40), label = c("10%","20%","30%", "40%") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  #geom_bar(aes(x=as.factor(id), y=perc, fill=disease), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  #geom_text(data=label_data, aes(x=id, y=perc+10, label=level_order, hjust=hjust), color="black",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=disease), hjust=c(1,1,1,1,0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p_genes

# png("/onek1k/Figures_and_Tables/Figure5A_overlap_genes_black.png",  units="px", width=2000, height=2000, res=300)
# p_genes
# dev.off()

# pdf("/onek1k/Figures_and_Tables/Figure5A_overlap_genes_white.pdf", height=8, width=8)
# p_genes
# dev.off()

pdf("/onek1k/Figures_and_Tables/Figure5A_overlap_genes_white_withoutlabels.pdf", height=8, width=8)
p_genes
dev.off()
