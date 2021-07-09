##############################################################################
# Script information                                                      
# Title: Conditioning plots per celltype pairs
# Author: Seyhan Yazar
# Date: 2021-03-15
# Description: None
##############################################################################

# import libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(car)
library(ggthemes)

workdir = '/onek1k/independent_lead_eQTL_analysis/cs_all_results'
setwd(workdir)

# get the filenames
filenames <- list.files(path=".", pattern="*_ALL_correlation_summary.tsv")

# read all the files
df_correlation <- lapply(filenames, function (x) { tryCatch(fread(x), error=function(e) NULL)})
length(df_correlation)
# df_correlation2 <- lapply (df_correlation, function(x) x <- x %>% group_by(geneid) %>% distinct)
name_lst <- sub("*_ALL_correlation_summary.tsv","",filenames)
names(df_correlation) <- name_lst
df_correlation2 <- bind_rows(df_correlation, .id="cells_to_compare")
colnames(df_correlation2)

# identify cell type 1 & 2 
df_correlation2$firstC <- sub("\\_.*","", df_correlation2$cells_to_compare)
df_correlation2$secondC <- sub(".*\\_","", df_correlation2$cells_to_compare)

cell_types_in_file <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown","NKmat", "NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC", "DC")
# cell_labels <- c("CD4 Naïve&CM","CD4 EM&TEMRA","CD4 SOX4","CD8 EM&TEMRA","CD8 Naïve&CM","CD8 S100B","NK", "NK Recruiting",
#                 "Plasma","B Mem","B Imm&Naïve","Mono C","Mono NC", "DC")

cell_labels <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET", "CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma", "B Mem", "B IN", "Mono C", "Mono NC", "DC")

cell_lst <- c(CD4all="CD4 NC", CD4effCM="CD4 ET", CD4TGFbStim="CD4 SOX4" ,
    CD8eff="CD8 ET", CD8all="CD8 NC", CD8unknown="CD8 S100B", 
    NKmat="NK", NKact="NK R", Plasma="Plasma" , 
    Bmem="B Mem", BimmNaive="B IN", MonoC="Mono C", MonoNC="Mono NC", DC="DC")

df_correlation2$celltype1 <- as.character(cell_lst[df_correlation2$firstC])
df_correlation2$celltype2 <- as.character(cell_lst[df_correlation2$secondC])

# identify significant snps after conditioning
df_correlation2$significance_status <- ifelse(df_correlation2$c1_new_p.value < 0.05, 1, 0)
df_correlation2$significance_status <- factor(df_correlation2$significance_status, levels=c(1,0), labels = c("significant", "non-significant"))

#### Create brand palette

# Color palette to use
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")

# Colour palette for cell types
branded_colors <- as.list(setNames(tol14rainbow,cell_labels))

  # Create a palette function
  branded_pal <- function(
  primary = "CD4 NC", 
  other = "grey", 
  direction = 1
) {
  stopifnot(primary %in% names(branded_colors))
  
  function(n) {
    if (n > 14 ) warning("Branded Color Palette only has 14 colors.")
    
    if (n == 2) {
      other <- if (!other %in% names(branded_colors)) {
        other
      } else {
        branded_colors[other]
      }
      color_list <- c(other, branded_colors[primary])
    } else {
      color_list <- branded_colors[1:n]
    }
    
    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}

# Create a discrete color scale
scale_colour_branded <- function(
  primary = "CD4 Naïve&CM", 
  other = "grey", 
  direction = 1, 
  ...
) {
  ggplot2::discrete_scale(
    "colour", "branded", 
    branded_pal(primary, other, direction), 
    ...
  )
}

scale_color_branded <- scale_colour_branded

dir.create("/onek1k/independent_lead_eQTL_analysis/plots")

####  Create graphing function ####
condition.graph <- function(df, na.rm = TRUE, ...){
  
  # create list of cell types in data to loop over 
  celltype_list <- unique(df$celltype1)
  celltype2_list <- unique(df$celltype2)

  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(celltype_list)) { 

    # create plot for each county in df 
    plot <- 
      ggplot(subset(df, df$celltype1==celltype_list[i]),
             aes(celltype1_estimate, c1_new_rho, group = celltype1, colour = significance_status)) + 
      geom_point(size=0.75) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      scale_colour_branded(primary=print(celltype_list[i]), direction= -1) +
      
      theme_pander(base_size=10, nomargin=FALSE) +
      theme(legend.position="none") +
      # , strip.background = element_blank(), strip.text.x = element_blank()) + 
      
      scale_y_continuous(paste('Correlation Coefficient of', celltype_list[i], 'Lead SNP After Regressing Lead eSNP of another cell type'), limits=c(-1,1)) +
      scale_x_continuous(paste("Original Correlation Coefficient (rho) of", celltype_list[i], "Lead SNP"),limits=c(-1,1)) +
      
      # facet_wrap( ~  celltype2, ncol=4) 
      facet_wrap( ~  celltype2, ncol=4, strip.position = "top", scales="free") 

      # ggtitle(paste(celltype_list[i])) 
    
    # save plots as .png
    ggsave(plot, file=paste('/onek1k/independent_lead_eQTL_analysis/conditioning_plots/',
                            celltype_list[i], ".png", sep=''), dpi="print")
   
    # save plots as .pdf
    #ggsave(plot, file=paste('/onek1k/independent_lead_eQTL_analysis/plots/',
    #                        celltype_list[i], ".pdf", sep=''), scale=2)
    
    # print plots to screen
    # print(plot)
  }
} 

condition.graph(df_correlation2)

# condition.graph(df_correlation2 %>% filter (celltype1=="CD8 Eff") %>% filter(celltype2=="Plasma"))

# plots <- lapply (dataset, function (id) {  
#    ggplot(data=, dataset$id, aes(x= id$celltype1_estimate, y=id$c1_new_rho)) +
#    geom_point() + 
#    theme_bw() +
#    scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(-1,1)) +
#    geom_hline(yintercept = 0, color = "black") +
#    geom_vline(xintercept = 0, color = "black")
# })
# m <- matrix(NA, 14, 14)
# m[lower.tri(m, diag = T)] <- 1:98
# grid_plots <- grid.arrange(grobs = plots, ncol=4)

# png("first_vs_second.png",units="px",res=300)
# grid_plots
# dev.off()

# ggsave("first_vs_second.png", grid_plots, dpi="print", width=800, height=800, units="mm")

df_for_corrplot <- df_correlation2 %>% 
  group_by(celltype1, celltype2, significance_status) %>% 
  summarise (n=n()) %>% 
  group_by(celltype1,celltype2) %>% 
  mutate(freq=n/sum(n)) %>% 
  filter(significance_status=="significant") %>%
  select(-significance_status, -n)

df_for_corrplot$corr <- cut(df_for_corrplot$freq, c(0, seq(0.1, 1, by=0.1)), include.lowest = T)

pdf(file="/onek1k/independent_lead_eQTL_analysis/plots/percentage_correlation_matrix.pdf")
df_for_corrplot %>%
  ggplot() + aes(celltype1, celltype2, fill = corr) + 
  geom_tile()
dev.off()





cell_labels <- c("CD4 Naïve&CM","CD4 EM&TEMRA","CD4 SOX4","CD8 EM&TEMRA","CD8 Naïve&CM","CD8 S100B","NK", "NK Recruiting",
                 "Plasma","B Mem","B Imm&Naïve","Mono C","Mono NC", "DC")

library(tidyr)

m_corrplot <- pivot_wider(df_for_corrplot, names_from="celltype1", values_from="freq") %>%  tibble::column_to_rownames("celltype2") %>% as.data.frame()
m_corrplot[is.na(m_corrplot)] = 0
m_corrplot <- m_corrplot[,cell_labels]
m_corrplot <- m_corrplot[cell_labels,]
m_corrplot <- as.matrix(m_corrplot)

library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

cols <- colorRamp2(c(0, 0.0001, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  c("white",brewer.pal(n=10, name="PuOr")))

pdf(file="/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/conditioning_analysis/plots/percentage_correlation_matrix.pdf")
ht <- Heatmap(m_corrplot, col = cols, 
  cluster_rows=FALSE, cluster_columns = FALSE,
  row_names_side = "left", column_names_side = "top",
  heatmap_legend_param = list(title="", at = c(0, 1), legend_height = unit(15, "cm"),  y = unit(3, "cm")))
ht 
dev.off()

### Additional Plots for Figure 3
# CD4 all vs B imm&Naive
df_example1 <- df_correlation2 %>% filter(celltype2=="CD4 NC") %>% filter(celltype1=="B IN")
table(df_example1$significance_status)
plot_example1 <- ggplot(df_example1, aes(celltype1_estimate, c1_new_rho, colour = significance_status)) +
      geom_point(size=2) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme_pander(base_size=14, nomargin=FALSE) +
      theme(legend.position="none") +
      scale_y_continuous(paste('Correlation Coefficient of B IN given CD4 NC'), limits=c(-.8,.8)) +
      scale_x_continuous(paste("Correlation Coefficient of B IN"),limits=c(-1,1)) +
      scale_color_grey(end = 0.8, start = 0.2)
ggsave(plot_example1, file=paste('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/conditioning_analysis/plots/example1.pdf'), scale=1)

# Plasma vs Mono nc
df_example2 <- df_correlation2 %>% filter(celltype1=="Plasma") %>% filter(celltype2=="Mono NC")
table(df_example2$significance_status)
plot_example2 <- ggplot(df_example2, aes(celltype1_estimate, c1_new_rho, colour = significance_status)) +
      geom_point(size=2) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme_pander(base_size=14, nomargin=FALSE) +
      theme(legend.position="none") +
      scale_y_continuous(paste('Correlation Coefficient of Plasma given Mono NC'), limits=c(-.8,.8)) +
      scale_x_continuous(paste("Correlation Coefficient of Plasma"),limits=c(-1,1)) +
      scale_color_grey(end = 0.8, start = 0.2)
ggsave(plot_example2, file=paste('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/conditioning_analysis/plots/example2.pdf'), scale=1)

# Check the allelic effect of SNP with change in direction 
df_example2_sign <- df_example2 %>% filter(significance_status=="significant")
dim(df_example2_sign)
df_example2_sign <- df_example2_sign %>% filter(celltype1_estimate > 0)
dim(df_example2_sign)
df_example2_sign <- df_example2_sign %>% filter(c1_new_rho < 0)
dim(df_example2_sign)


# Mono nc vs B IN
df_example2 <- df_correlation2 %>% filter(celltype1=="Mono NC") %>% filter(celltype2=="B IN")
table(df_example2$significance_status)
plot_example2 <- ggplot(df_example2, aes(celltype1_estimate, c1_new_rho, colour = significance_status)) +
      geom_point(size=2) + 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme_pander(base_size=14, nomargin=FALSE) +
      theme(legend.position="none") +
      scale_y_continuous(paste('Correlation Coefficient of Mono NC given B IN'), limits=c(-.8,.8)) +
      scale_x_continuous(paste("Correlation Coefficient of Mono NC"),limits=c(-1,1)) +
      scale_color_grey(end = 0.8, start = 0.2)
ggsave(plot_example2, file=paste('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/conditioning_analysis/plots/example3.pdf'), scale=1)
