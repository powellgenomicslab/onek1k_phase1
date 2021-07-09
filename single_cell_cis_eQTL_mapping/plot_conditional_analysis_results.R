##############################################################################
# Script information                                                      
# Title: Conditional Analysis Results - Figure 2
# Author: Seyhan Yazar
# Date: 2020-12-25
# Description: This R script was written to generate Figure 2 plots
##############################################################################

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

# Add these rsIDs 
dataset2[dataset2$snpid=="17:44355634_A", "ID"] <- "rs2696531"

# Function to select all esnps for each cell type
nth_element <- function(vector, starting_position) {
  vector[seq(starting_position,starting_position+4)]
}

# Combine results by cell type
create_df_by_celltype <- function (first.file.number) {
  cellFiles = nth_element (files,first.file.number)
  print(cellFiles)
  cellFiles_df = lapply(cellFiles, function(x) tryCatch(read_tsv(x), error=function(e) NULL)) %>% bind_rows() 
  cellFiles_df
}
first_file_numbers <- seq(1, length(files),5)
dfList_by_celltype <- lapply(first_file_numbers, create_df_by_celltype)
firstFilenames = files[first_file_numbers]
celltypes = sub("*_eSNP1.tsv", "", firstFilenames)
names(dfList_by_celltype) <- celltypes

dataset2$celltype <- factor(dataset2$celltype, levels=cell_order, 
                            labels=corder)

# # identify the list of SNPs
snp <- unlist(lapply (dfList_by_celltype, function(x) x%>% select(snpid)),recursive = F)
snp2 <-data.frame(unique(as.vector(unlist(snp, recursive = F, use.names = F))))
colnames(snp2)[1] <-"snpid"
snp2$snps <- as.character(snp2$snpid)

# Create the upset plot dataset
dataset_upset <- data.frame(do.call(cbind,lapply (dfList_by_celltype, function (x) ifelse(snp2$snpid %in% x$snpid,1,0))))
dataset_upset %>% head()
dataset_upset2 <- dataset_upset[rev(cell_order)]
dataset_upset2 %>% head()
labels=rev(corder)
colnames(dataset_upset2) <- labels
dataset_upset2 %>% head()

# Rankings - Figure 2A
dataset2$celltype <- factor(dataset2$celltype, levels=corder, labels= corder)
pdf("/onek1k/Figures_and_Tables/Figure2A_conditional_eQTLs_by_rank.pdf", width=9, height=9)
dataset2 %>% 
  group_by(corder, eSNP_rank) %>% count(eSNP_rank)  %>%
  mutate(corder=factor(corder, levels=cells$corder)) %>%
  ggplot(aes(x=eSNP_rank, y=n, fill=corder)) + 
  geom_bar(stat="identity") +
  facet_wrap(.~corder, ncol=4, scales = "free") +
  theme_classic(15) +
  theme(legend.position = "none") +
  scale_fill_manual(values = tol14rainbow) +
  theme(strip.background = element_blank()) +
  ylab("cis-eQTLs") +
  scale_x_discrete(name ="cis-eQTL ranks", 
                   labels=c("1","2","3","4","5")) 
dev.off()  

# Upset Plot Simple - Figure 2B
pdf("/onek1k/Figures_and_Tables/Figure2B_UpsetPlot.pdf",  
    width=15, height=8)
upset(dataset_upset2,  sets = colnames(dataset_upset2), keep.order = TRUE,
      matrix.color = "darkblue", main.bar.color = "gray5", mainbar.y.label = "Independent eQTLs",
      sets.bar.color = rev(tol14rainbow), sets.x.label = "Total eQTLs (10^3)", point.size = 3,
      mb.ratio = c(0.6, 0.4), line.size = 0.3, order.by = c("freq","degree"), text.scale = 2,
      show.numbers ="no") 
dev.off()

# Identify unique genes
unique_genes<- unlist(lapply (dfList_by_celltype, function(x) x%>% select(geneid)),recursive = F)
unique_genes2 <-data.frame(unique(as.vector(unlist(unique_genes, recursive = F, use.names = F))))
colnames(unique_genes2)[1] <-"geneid"
unique_genes2$geneid <- as.character(unique_genes2$geneid)


### Expression of eGenes - Figure 2C
colnames(unique_genes2)[1] <-"gene"
expression_levels <- read.table("/onek1k/perc_means_all.tsv",header = T,se="\t")
colnames(expression_levels)[2] <- "cell_labels"
expression_levels <- expression_levels[(expression_levels$gene %in% unique_genes2$gene),]
expression_levels$gene <- factor(expression_levels$gene)
df_exp <- expression_levels[c("gene", "cell_labels","means")]
df_exp <- df_exp[df_exp$cell_labels %in% c("CD4 all", "CD4 Eff & CM","CD4 TGFb stim",
                                           "CD8 all", "CD8 Eff", "CD8 unknown",
                                           "NK mat", "NK act", "Plasma",
                                           "B Mem", "B Imm & Naive",
                                           "Mono c", "Mono nc", "DC"),]
sub_df_exp <- df_exp %>% filter(cell_labels %in% c("Mono c", "Plasma"))
nrow(sub_df_exp)
sub_df_exp %>% count(gene)

m_exp <- pivot_wider(df_exp, names_from = "cell_labels",values_from = "means") %>% column_to_rownames("gene") %>% as.data.frame()
m_exp[is.na(m_exp)] = 0
colnames (m_exp) <- cells$corder
m_exp <- m_exp[,cell_labels]
m_exp <- m_exp[,corder]

cpairs <- expand_grid(i=colnames(m_exp), j=colnames(m_exp))
cpairs2 <- subset(cpairs, i!=j)

find_concordance <- function (x) {
  c1 <- cpairs2$i[x]
  c2 <- cpairs2$j[x]
  data <- m_exp %>% select(print(c1), print(c2))
  coefficent <- epi.ccc(data[,1],data[,2])$rho.c[,1]
  output <- data.frame(c1,c2,coefficent)
  output
}

expr_coeffs <- lapply(1:182, find_concordance)
expr_coeffs_df <- do.call(rbind.data.frame, expr_coeffs)
m_corr <- pivot_wider(expr_coeffs_df, names_from = "c1",values_from = "coefficent") %>% column_to_rownames("c2") %>% as.data.frame()
m_corr  <- m_corr[,rev(corder)]
m_corr  <- m_corr[corder,]
m_corr <- as.matrix(m_corr)
m_corr[is.na(m_corr)] = 1
p <- ggcorrplot::ggcorrplot(m_corr)
pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance.pdf",  width=12, height=9)
p + scale_fill_gradient2(limit = c(0.75,1), low = "#B2ABD2", high = "#2D004B", mid = "#FEE0B6", midpoint = 0.875,
                         name="Lin's concordance\ncorrelation coefficient") + theme(legend.position = "right") + 
  guides(fill=guide_legend(label.theme = element_text(angle = 360))) 
dev.off()

cols <- colorRamp2(c(0.75, 0.8, 0.85, 0.9, 0.95, 1 ),  c("white",brewer.pal(n=5, name="PuOr")))

ht <- Heatmap(m_corr, col = cols, 
              cluster_rows=FALSE, cluster_columns = FALSE,
              column_order = corder,
              row_names_side = "left", column_names_side = "top",
              heatmap_legend_param = list(title="Lin's concordance \ncorrelation coefficietn", 
                                          legend_height = unit(15, "cm"),  y = unit(3, "cm"), 
                                          at = c(0.75, 0.80, 0.85, 0.90, 0.95, 1)))
ht 

pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance_210409.pdf",  width=11, height=9)
ht
dev.off()

# Example 1 
m_exp_exp <- m_exp %>% select("CD4 NC","CD8 NC")
tmp.ccc <- epi.ccc(m_exp_exp[,1], m_exp_exp[,2], ci = "z-transform", 
                   conf.level = 0.95, rep.measure = FALSE)
tmp.lab <- data.frame(lab = paste("CCC: ", 
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ", 
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))
z <- lm(m_exp_exp[,1] ~ m_exp_exp[,2])
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance_example1.pdf",  
    width=9, height=9)
ggplot(m_exp_exp, aes(x = m_exp_exp[,1], y = m_exp_exp[,2])) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta), colour = "#2D004B", size = 2) +
  # xlim(0, 5) +
  # ylim(0, 5) +
  xlab("mean expression in CD4 NC") +
  ylab("mean expression in CD8 NC") +
  geom_text(data = tmp.lab, x = 1.5, y = 4.9, label = tmp.lab$lab, fontface = "bold", colour = "#2D004B") + 
  coord_fixed(ratio = 1 / 1) +
  theme_minimal(30)
dev.off()

# Example 2

m_exp_exp <- m_exp %>% select('Mono NC','Plasma')
tmp.ccc <- epi.ccc(m_exp_exp[,1], m_exp_exp[,2], ci = "z-transform", 
                   conf.level = 0.95, rep.measure = FALSE)
tmp.lab <- data.frame(lab = paste("CCC: ", 
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ", 
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))
z <- lm(m_exp_exp[,1] ~ m_exp_exp[,2])
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance_example.pdf",  
    width=9, height=9)
ggplot(m_exp_exp, aes(x = m_exp_exp[,1], y = m_exp_exp[,2])) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta), colour = "#7F3B08", size = 2) +
  xlim(0, 5) +
  ylim(0, 5) +
  xlab("mean expression in Mono NC") +
  ylab("mean expression in Plasma") +
  geom_text(data = tmp.lab, x = 1.5, y = 4.9, label = tmp.lab$lab, fontface = "bold", colour = "#7F3B08") + 
  coord_fixed(ratio = 1 / 1) +
  theme_minimal(30)
dev.off()

# Example 3

m_exp_exp <- m_exp %>% select('Mono NC','B IN')
tmp.ccc <- epi.ccc(m_exp_exp[,1], m_exp_exp[,2], ci = "z-transform", 
                   conf.level = 0.95, rep.measure = FALSE)
tmp.lab <- data.frame(lab = paste("CCC: ", 
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ", 
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))
z <- lm(m_exp_exp[,1] ~ m_exp_exp[,2])
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance_example3.pdf",  
    width=9, height=9)
ggplot(m_exp_exp, aes(x = m_exp_exp[,1], y = m_exp_exp[,2])) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta), colour = "#7F3B08", size = 2) +
  xlim(0, 5) +
  ylim(0, 5) +
  xlab("mean expression in Mono NC") +
  ylab("mean expression in B IN") +
  geom_text(data = tmp.lab, x = 1.5, y = 4.9, label = tmp.lab$lab, fontface = "bold", colour = "#7F3B08") + 
  coord_fixed(ratio = 1 / 1) +
  theme_minimal(30)
dev.off()

# Example 3

m_exp_exp <- m_exp %>% select('Mono C','Plasma')
tmp.ccc <- epi.ccc(m_exp_exp[,1], m_exp_exp[,2], ci = "z-transform", 
                   conf.level = 0.95, rep.measure = FALSE)
tmp.lab <- data.frame(lab = paste("CCC: ", 
                                  round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ", 
                                  round(tmp.ccc$rho.c[,2], digits = 2), " - ",
                                  round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))
z <- lm(m_exp_exp[,1] ~ m_exp_exp[,2])
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

pdf("/onek1k/Figures_and_Tables/Figure2C_exprs_concordance_example4.pdf",  
    width=9, height=9)
ggplot(m_exp_exp, aes(x = m_exp_exp[,1], y = m_exp_exp[,2])) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta), colour = "#7F3B08", size = 2) +
  xlim(0, 5) +
  ylim(0, 5) +
  xlab("mean expression in Mono C") +
  ylab("mean expression in Plasma") +
  geom_text(data = tmp.lab, x = 1.5, y = 4.9, label = tmp.lab$lab, fontface = "bold", colour = "#7F3B08") + 
  coord_fixed(ratio = 1 / 1) +
  theme_minimal(30)
dev.off()
