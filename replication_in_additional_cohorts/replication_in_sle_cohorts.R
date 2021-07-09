##############################################################################
# Script information                                                      
# Title: Replication in SLE cohorts
# Author: Seyhan Yazar
# Date: 2021-03-12
# Description: Includes figure 3 plots
##############################################################################

# Set output ------------------------------------------
output <- set_output("2021-03-12", "replication-in-sle")

# Import libraroes
library("R.utils")
library("data.table")
library("dplyr")
library("tidyr")
library("ComplexHeatmap")
library("circlize")
library("ggplot2")
library("shades")
library(stringr)
library("RColorBrewer")

#####
# Colour palettee
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
# Order of cells
celltype <- c("CD4 NC","CD4 ET","CD4 SOX4","CD8 ET","CD8 NC","CD8 S100B","NK", "NK R",
                 "Plasma","B Mem","B IN","Mono C","Mono NC", "DC")

cell_type <- c("CD4all", "CD4effCM", "CD4TGFbStim", "CD8eff", "CD8all", "CD8unknown", "NKmat", "NKact", "Plasma", "Bmem",
              "BimmNaive", "MonoC", "MonoNC", "DC")

cells <- data.frame(cell_type, celltype)

#### OneK1K Study ####
onek1k_df <- fread("/onek1k/OneK1K_eQTLs_Results_03022021.tsv")
head(onek1k_df)
table(onek1k_df$cell_type)
onek1k_df$eSNP <- paste0(onek1k_df$Chromosome,":", onek1k_df$Position, "_", onek1k_df$SNP_assessed_allele)
onek1k_df <- onek1k_df %>% filter(eSNP_rank=="eSNP1")

# onek1k_df <- onek1k_df %>% filter(cell_type!="Plasma")
colnames(onek1k_df)
onek1k_df <- onek1k_df[,c("cell_type", "GeneID", "eSNP", "rho_correlation_coefficient")]
colnames(onek1k_df) <- c("cell_type", "geneid", "snpid", "onek1k_estimate")

onek1k_df_b <- onek1k_df %>% filter(cell_type %in% c("B IN", "B Mem", "Plasma"))
onek1k_df_b_p <- onek1k_df_b %>% select("snpid","geneid", "onek1k_estimate")

onek1k_df_cd4 <- onek1k_df %>% filter(cell_type %in% c("CD4 NC", "CD4 ET", "CD4 SOX4"))
onek1k_df_cd4_p <- onek1k_df_cd4 %>% select("snpid","geneid", "onek1k_estimate")

onek1k_df_cd8 <- onek1k_df %>% filter(cell_type %in% c("CD8 ET", "CD8 NC", "CD8 S100B"))
onek1k_df_cd8_p <- onek1k_df_cd8 %>% select("snpid","geneid", "onek1k_estimate")

onek1k_df_dc <- onek1k_df %>% filter(cell_type %in% c("DC"))
onek1k_df_dc_p <- onek1k_df_dc %>% select("snpid","geneid", "onek1k_estimate")

onek1k_df_mono <- onek1k_df %>% filter(cell_type %in% c("Mono C", "Mono NC"))
onek1k_df_mono_p <- onek1k_df_mono %>% select("snpid","geneid", "onek1k_estimate")

onek1k_df_nk <- onek1k_df %>% filter(cell_type %in% c("NK", "NK R"))
onek1k_df_nk_p <- onek1k_df_nk %>% select("snpid","geneid", "onek1k_estimate")

###### Replication Cohorts ######
setwd("/onek1k/correlation_outputs_210423")
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
sle_df <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(sle_df) <- filenames
sle_df <- bind_rows(sle_df, .id="filename")
sle_df$oldname <- gsub("_.*","",sle_df$filename)
sle_df$ye_celltype <- str_match(sle_df$filename, "eSNP1_\\s*(.*?)\\s*_correlation")[,2]
sle_df$cohort <- factor(sub(".*_", "", sle_df$filename))

# Regardless significance - all tested
ye_rep_df <- sle_df %>% filter(cohort!=c("immvar"))
ye_rep_df$significant <- ifelse(ye_rep_df$localFDR < 0.1, "yes", "no")
colnames(ye_rep_df)[8] <- "cell_type"

ye_df_b <- ye_rep_df %>% filter(cell_type %in% c("BimmNaive", "Bmem"))
ye_df_b_p <- ye_df_b %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_cd4 <- ye_rep_df %>% filter(cell_type %in% c("CD4all", "CD4effCM", "CD4TGFbStim"))
ye_df_cd4_p <- ye_df_cd4 %>% select("snpid","geneid", "estimate","significant", "cohort" )

ye_df_cd8 <- ye_rep_df %>% filter(cell_type %in% c("CD8all", "CD8eff", "CD8unknown"))
ye_df_cd8_p <- ye_df_cd8 %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_dc <- ye_rep_df %>% filter(cell_type %in% c("DC"))
ye_df_dc_p <- ye_df_dc %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_mono <-ye_rep_df %>% filter(cell_type %in% c("MonoC", "MonoNC"))
ye_df_mono_p <- ye_df_mono %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_nk <- ye_rep_df %>% filter(cell_type %in% c("NKact", "NKmat"))
ye_df_mono_p  <- ye_df_nk %>% select("snpid","geneid", "estimate", "significant", "cohort")

# Overlap function
overlap <- function(df, one) {
  rep_df <- left_join(df, one, by=c("snpid", "geneid"))
  rep_df <- rep_df[complete.cases(rep_df), ]
  rep_df <- rep_df %>% mutate (onek1k_direct = ifelse(onek1k_estimate < 0 , 0 , 1 ))
  rep_df <- rep_df %>% mutate (replication = ifelse(estimate < 0, 0 , 1 ))
  rep_df$agreement <- rep_df$onek1k_direct + rep_df$replication
  rep_df$agreement
  rep_df$agreement <- factor(rep_df$agreement, levels=c("0", "1", "2"), labels = c("yes", "no", "yes"))
  rep_df
}

rep_df_b_v1 <- overlap(ye_df_b_p, onek1k_df_b_p)
rep_df_cd4_v1 <- overlap(ye_df_cd4_p, onek1k_df_cd4_p)
rep_df_cd8_v1 <- overlap(ye_df_cd8_p, onek1k_df_cd8_p)
rep_df_dc_v1 <- overlap(ye_df_dc_p , onek1k_df_dc_p)
rep_df_mono_v1 <- overlap(ye_df_mono_p ,onek1k_df_mono_p)
rep_df_nk_v1  <- overlap(ye_df_mono_p , onek1k_df_nk_p)
results_list <- list(rep_df_b_v1, rep_df_cd4_v1, rep_df_cd8_v1, rep_df_dc_v1,
                     rep_df_mono_v1, rep_df_nk_v1)
names(results_list) <- c("B-cell", "CD4-T-cell","CD8-T-cell",
                         "DC", "Monocyte", "NK")

combined_df <- bind_rows(results_list, .id="cell_type")

new_df <- combined_df %>% select(cell_type, agreement, cohort)
new_df2 <- new_df %>% group_by(cohort, cell_type, agreement) %>% summarise(n = n()) %>%
  mutate(freq = (n / sum(n))*100)
new_df2$agreement <- factor(new_df2$agreement, levels=rev(c("yes", "no")),labels=c("Disconcordant", "Concordant"))
new_df2$cell_type <- factor(new_df2$cell_type, levels=c("CD4-T-cell", "CD8-T-cell",
                                                        "NK", "B-cell","Monocyte", "DC"))

brewer.pal(12, "Paired")
direction_plot <- new_df2 %>%
  ggplot(aes(y=freq, x=cohort, fill=agreement)) +
  geom_bar(stat="identity", position="stack") +
  geom_hline(yintercept = 50, linetype="dashed") +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values=c("#FDB863","#B35806")) +
  ggthemes::theme_clean(15) +
  theme(axis.title = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse=TRUE)) +
  facet_grid(cell_type ~.)
pdf("/onek1k/Figures_and_Tables/Figure3A_Replication_all_tested.pdf",
    width=7, height=10)
direction_plot
dev.off()

# Replication @5% FDR
ye_rep_df <- sle_df %>% filter(cohort!=c("immvar"))
ye_rep_df$significant <- ifelse(ye_rep_df$localFDR < 0.1, "yes", "no")
ye_rep_sig_df <-ye_rep_df %>% filter(significant=="yes")
colnames(ye_rep_sig_df)[8] <- "cell_type"

ye_df_b_sig <- ye_rep_sig_df %>% filter(cell_type %in% c("BimmNaive", "Bmem"))
ye_df_b_p_sig <- ye_df_b_sig %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_cd4_sig <- ye_rep_sig_df %>% filter(cell_type %in% c("CD4all", "CD4effCM", "CD4TGFbStim"))
ye_df_cd4_p_sig <- ye_df_cd4_sig %>% select("snpid","geneid", "estimate","significant", "cohort" )

ye_df_cd8_sig <- ye_rep_sig_df %>% filter(cell_type %in% c("CD8all", "CD8eff", "CD8unknown"))
ye_df_cd8_p_sig <- ye_df_cd8_sig %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_dc_sig <- ye_rep_sig_df %>% filter(cell_type %in% c("DC"))
ye_df_dc_p_sig <- ye_df_dc_sig %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_mono_sig <-ye_rep_sig_df %>% filter(cell_type %in% c("MonoC", "MonoNC"))
ye_df_mono_p_sig <- ye_df_mono_sig %>% select("snpid","geneid", "estimate", "significant", "cohort")

ye_df_nk_sig <- ye_rep_sig_df %>% filter(cell_type %in% c("NKact", "NKmat"))
ye_df_nk_p_sig  <- ye_df_nk_sig %>% select("snpid","geneid", "estimate", "significant", "cohort")

# Overlap function
rep_df_b <- overlap(ye_df_b_p_sig, onek1k_df_b_p)
rep_df_cd4 <- overlap(ye_df_cd4_p_sig, onek1k_df_cd4_p)
rep_df_cd8 <- overlap(ye_df_cd8_p_sig, onek1k_df_cd8_p)
rep_df_dc <- overlap(ye_df_dc_p_sig , onek1k_df_dc_p)
rep_df_mono <- overlap(ye_df_mono_p_sig ,onek1k_df_mono_p)
rep_df_nk  <- overlap(ye_df_mono_p_sig , onek1k_df_nk_p)
results_list <- list(rep_df_b, rep_df_cd4, rep_df_cd8, rep_df_dc,
                     rep_df_mono, rep_df_nk)
names(results_list) <- c("B-cell", "CD4-T-cell","CD8-T-cell",
                         "DC", "Monocyte", "NK")
results_df <- bind_rows(results_list, .id="cell_type")

new_df <- results_df %>% select(cell_type, agreement, cohort)
new_df2 <- new_df %>% group_by(cohort, cell_type, agreement) %>% summarise(n = n()) %>%
  mutate(freq = (n / sum(n))*100)
new_df2$agreement <- factor(new_df2$agreement, levels=rev(c("yes", "no")),labels=c("Disconcordant", "Concordant"))
new_df2$cell_type <- factor(new_df2$cell_type, levels=c("CD4-T-cell", "CD8-T-cell",
                                                            "NK", "B-cell","Monocyte", "DC"))

brewer.pal(12, "Paired")
direction_plot <- new_df2 %>%
  ggplot(aes(y=freq, x=cohort, fill=agreement)) +
  geom_bar(stat="identity", position="stack") +
  geom_hline(yintercept = 50, linetype="dashed") +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values=c("#FDB863","#B35806")) +
  ggthemes::theme_clean(15) +
  theme(axis.title = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank()) +
  guides(fill = guide_legend(reverse=TRUE)) +
  facet_grid(cell_type ~.)
pdf("/onek1k/Figures_and_Tables/Figure3B_Replication_5percFDR.pdf",
    width=7, height=10)
direction_plot
dev.off()

### Figure 2C Overlap of eGenes between cell types
eQTL_euro<- sle_df %>% filter(cohort==c("euro"))
index <- which(eQTL_euro$localFDR<0.05)
genes_euro <- as.data.frame(table(table(eQTL_euro$geneid[index])))
names(genes_euro)[1] <- "CellTypes"
names(genes_euro)[2] <- "Genes"
genes_euro$Cohort <- "EUR"

eQTL_asian<- sle_df %>% filter(cohort==c("asian"))
index <- which(eQTL_asian$localFDR<0.05)
genes_asian <- as.data.frame(table(table(eQTL_asian$geneid[index])))
names(genes_asian)[1] <- "CellTypes"
names(genes_asian)[2] <- "Genes"
genes_asian$Cohort <- "EAS"

genes <- rbind(genes_euro, genes_asian)
genes$Cohort <- factor(genes$Cohort, levels=c("EUR", "EAS"))

pdf("/onek1k/Figures_and_Tables/Figure3C_Replication_eGenes.pdf",
    height=5, width=8)
ggplot(genes) + geom_bar(stat="identity", aes(x=CellTypes, y=Genes, fill=Cohort), position=position_dodge()) + 
  #coord_flip() +
  theme_classic(15) +
  xlab("Number of cell types") +
  ylab("Number of eGenes") +
  scale_fill_manual(values=c("orange", "darkorange3"))
dev.off()

### Figure 2D Overlap of eSNPs between cell types
eSNPs_euro <- as.data.frame(table(table(eQTL_euro$snpid[index])))
names(eSNPs_euro)[1] <- "CellTypes"
names(eSNPs_euro)[2] <- "eSNPs"
eSNPs_euro$Cohort <- "EUR"

eSNPs_asian <- as.data.frame(table(table(eQTL_asian$snpid[index])))
names(eSNPs_asian)[1] <- "CellTypes"
names(eSNPs_asian)[2] <- "eSNPs"
eSNPs_asian$Cohort <- "EAS"

eSNPs <- rbind(eSNPs_euro, eSNPs_asian)
eSNPs$Cohort <- factor(eSNPs$Cohort, levels=c("EUR", "EAS"))

pdf("/onek1k/Figures_and_Tables/Figure3D_Replication_eSNPs.pdf",
    height=5, width=8)
ggplot(eSNPs) + geom_bar(stat="identity", aes(x=CellTypes, y=eSNPs, fill=Cohort), position=position_dodge()) +  
  #coord_flip() +
  theme_classic(15) +
  xlab("Number of cell types") +
  ylab("Number of eSNPs") +
  scale_fill_manual(values=c("#b53c6d","#892D52"))
dev.off()

####
combined_df$cell_type <- factor(combined_df$cell_type, levels=c("CD4-T-cell", "CD8-T-cell",
                                                            "NK", "B-cell","Monocyte", "DC"))
pdf("/onek1k/Figures_and_Tables/Figure3E_European.pdf",
    width=14, height=3)
combined_df %>% 
  filter (cohort=="euro") %>%
  group_by(cell_type) %>%
  mutate(significant = if_else(as.character(significant) == "yes", as.character(cell_type), "zzz")) %>% 
  mutate(significant = as.factor(significant)) %>% 
  ggplot() +
  aes(x=onek1k_estimate, y=estimate, color=significant) +
  geom_point(size=1) +
  facet_wrap(~cell_type, ncol=6) +
  theme_minimal(15) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  scale_colour_manual(values=c("#882E72", "grey","#1965B0", "#4EB265", "#F7EE55","#F1932D", "#DC050C")) +
  # scale_color_manual(values = brewer.pal(n=6, name="Set2")) +
  xlab("Correlation coefficient in OneK1K Study") +
  ylab("Correlation coefficient \nin Replication Cohort  \n(European samples)")
dev.off()

pdf("/onek1k/Figures_and_Tables/Figure3E_Asian.pdf",
    width=14, height=3)
combined_df %>% 
  filter (cohort=="asian") %>%
  group_by(cell_type) %>%
  mutate(significant = if_else(as.character(significant) == "yes", as.character(cell_type), "zzz")) %>% 
  mutate(significant = as.factor(significant)) %>% 
  ggplot() +
  aes(x=onek1k_estimate, y=estimate, color=significant) +
  geom_point(size=1) +
  facet_wrap(~cell_type, ncol=6) +
  theme_minimal(15) +
  theme(legend.position = "none") +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  scale_colour_manual(values=c("#882E72", "grey","#1965B0", "#4EB265", "#F7EE55","#F1932D", "#DC050C")) +
  # scale_color_manual(values = brewer.pal(n=6, name="Set2")) +
  xlab("Correlation coefficient in OneK1K Study") +
  ylab("Correlation coefficient \nin Replication Cohort  \n(Asian samples)")
dev.off()

length(unique(onek1k_df_b$snpid))
length(unique(onek1k_df_cd4$snpid))
length(unique(onek1k_df_cd8$snpid))
length(unique(onek1k_df_dc$snpid))
length(unique(onek1k_df_mono$snpid))
length(unique(onek1k_df_nk$snpid))

ye_rep_df$cell_type <- ye_rep_df$oldname
ye_rep_df <- left_join(ye_rep_df, cells, by="cell_type")
ye_rep_df <- ye_rep_df %>% select("celltype", "geneid", "snpid", "estimate", "statistic", "p.value", "localFDR", "cohort")
ye_rep_df_euro <- ye_rep_df %>% filter(cohort=="euro") %>% select(-c("cohort"))
fwrite(ye_rep_df_euro, "/onek1k/Figures_and_Tables/TableSX_EUR.tsv", sep="\t")
ye_rep_df_asian <- ye_rep_df %>% filter(cohort=="asian") %>% select(-c("cohort"))
fwrite(ye_rep_df_asian, "/onek1k/Figures_and_Tables/TableSX_EAS.tsv", sep="\t")

# Session info ------------------------------------------
print_session(here(output))