# Script information ------------------------------------------------------

# title: Create plots for cell type proportions
# author: Jose Alquicira Hernandez
# date: 2020/10/01
# description: None

# qrsh -N prop -l mem_requested=200G,tmp_requested=200G 
# conda activate r-4.0.2 

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")

# Set output --------------------------------------------------------------

output <- set_output("2020-10-01", "cell_type_proportions")

# Read data ---------------------------------------------------------------

# Input
path <- file.path("results", "2019-08-13_cell_type_anno")
filename <- "cell_type.RDS"
input <-  here(path, filename)

# Read file
inicio("Reading gene expression data")
data <- readRDS(input)
fin()

# Remove vst output
data@assays$SCT@misc$vst.out <- NULL

# Remove biological outliers ----------------------------------------------
inicio("Detect biological outliers")
cells <- !data$cell_type %in% c("Platelets", "Erythrocytes")
data <- data[, cells]
fin()


data@meta.data %>% 
  as_tibble() %>% 
  group_by(cell_type, individual) %>% 
  tally() -> total


total %>% 
  pivot_wider(id_cols = "cell_type", names_from = "individual", 
              values_from = "n") %>% 
  as.data.frame() %>% 
  column_to_rownames("cell_type") -> total_wide

map_df(total_wide, ~{
  .x / sum(.x) * 100
}) %>% as.data.frame() -> total_wide_perc

rownames(total_wide_perc) <- rownames(total_wide)
stats <- apply(total_wide_perc, 1, function(x){
  x <- x[!is.na(x)]
  summary(x)
})
stats <- t(stats)
stats <- round(stats, 2)
stats <- as.data.frame(stats)
stats <- rownames_to_column(stats, "Cell type")
write_csv(stats, here(output, "props.csv"))


head(data[[]])


data[[]] %>%
  select(nCount_RNA, cell_type) %>% 
  group_by(cell_type) %>% 
  summarise_at(vars(nCount_RNA), list(min = min, 
                                      q25 = ~quantile(., 0.25), 
                                      median = median, 
                                      mean = mean,
                                      q75 = ~quantile(., 0.75), 
                                      max = max
  )
  ) -> ncounts

write_csv(ncounts, here(output, "counts.csv"))



data[[]] %>%
  select(nCount_SCT, cell_type) %>% 
  group_by(cell_type) %>% 
  summarise_at(vars(nCount_SCT), list(min = min, 
                                      q25 = ~quantile(., 0.25), 
                                      median = median, 
                                      mean = mean,
                                      q75 = ~quantile(., 0.75), 
                                      max = max
  )
  ) -> ncounts

write_csv(ncounts, here(output, "counts_sct.csv"))

# Relabel cells -----------------------------------------------------------

# Drop platelet and red cell levels if any
data$cell_type <- fct_drop(data$cell_type)
# Set order of plotting
level_order <- c("CD4 NC", "CD4 TE", "CD4 SOX4", 
                 "CD8 TE", "CD8 NC","CD8 S100B", 
                 "NK", "NK R", 
                 "Plasma", "B Mem", "B IN",  
                 "Mono C", "Mono NC", 
                 "DC")
# Map labels
cell_type <- recode(data$cell_type,
                    `CD4+ KLRB1+ T cell`      = "CD4 TE",
                    `CD4+ KLRB1- T cell`      = "CD4 NC",
                    `CD4+ SOX4+ T cell`       = "CD4 SOX4",
                    `CD8+ LTB+ T cell`        = "CD8 NC",
                    `CD8+ GNLY+ NKG7+ T cell` = "CD8 TE",
                    `CD8+ S100B+ T cell`      = "CD8 S100B",
                    `XCL1- NK`                = "NK",
                    `XCL1+ NK`                = "NK R",
                    `TCL1A- FCER2- B cell`    = "B Mem",
                    `TCL1A+ FCER2+ B cell`    = "B IN",
                    `IgJ+ B cell`             = "Plasma",
                    `Monocyte CD14+`          = "Mono C",
                    `Monocyte FCGR3A+`        = "Mono NC",
                    `Dendritic cell`          = "DC"
) 
# Create new column with new labels
data$new_cell_type <- factor(cell_type, level_order)


# Distribution of cell type proportions in 75 sequencing pools  -----------
pal <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", 
         "#E8601C", "#DC050C")


data[[]] %>% 
  as_tibble() %>% 
  mutate(pool = as.character(pool)) %>% 
  group_by(pool, new_cell_type) %>% 
  tally() %>% 
  ungroup %>% 
  group_by(pool) %>% 
  mutate(perc = n / sum(n) * 100) %>% 
  mutate(pool = str_remove(pool, "pool_")) %>% 
  ungroup() %>% 
  mutate(pool = factor(pool, pool %>% unique() %>% as.numeric() %>% sort() %>% as.character())) %>%
  ggplot() +
  aes(x = pool, y = perc, fill = new_cell_type) +
  geom_bar(stat = "identity", width = 1, color = "black", lwd = 0.1) +
  xlab("Pool") +
  ylab("Percentage of cells (%)") +
  theme_pub() +
  rotate_x() +
  guides(fill = guide_legend(title = "Cell type")) +
  scale_fill_manual(values = pal) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 18)) -> p_pools


ggsave(here(output, "cell-type_pool.png"), p_pools, width = 16, height = 6, dpi = "print")
ggsave(here(output, "cell-type_pool.pdf"), p_pools, width = 14, height = 6)


# Distribution of cell percentages among 982 individuals in 14 cell types -

data[[]] %>% 
  as_tibble() %>% 
  group_by(individual, new_cell_type) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(individual = str_split(individual, "_") %>% map_chr(2)) %>% 
  group_by(individual) %>% 
  mutate(perc = n / sum(n) * 100) %>% 
  ungroup() %>% 
  select(-n) %>% 
  as.data.frame() ->  prop_ind


prop_ind %>% 
  as_tibble() %>% 
  group_by(new_cell_type) %>% 
  mutate(individual = fct_reorder(individual, perc, mean)) -> prop_ind

ggplot(prop_ind) +
  aes(y = individual, x = perc, fill = new_cell_type) +
  geom_bar(stat = "identity", width = 1, color = "black", lwd = 0.05) +
  xlab("Percentage of cells (%)") +
  ylab("Individual") +
  theme_pub() +
  rotate_x() +
  guides(fill = guide_legend(title = "Cell type")) +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 18)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> p_indvs


ggsave(here(output, "cell-type_individual.png"), p_indvs, height = 16, width = 6, dpi = "print")
ggsave(here(output, "cell-type_individual.pdf"), p_indvs,height = 16, width = 6, dpi = "print")


# Distribution of cell percentages among 982 individuals in 14 cell types -----

library(ggridges)

prop_ind %>% 
  mutate(new_cell_type = fct_rev(new_cell_type)) %>% 
ggplot() +
  aes(x = perc, y = new_cell_type, fill = new_cell_type, color = new_cell_type) +
  geom_density_ridges(lwd = 0.05,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7) +
  xlab("Percentage of cells (%)") +
  ylab("Cell type") +
  theme_pub() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  theme(axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.position = "none") -> p_perc


ggsave(here(output, "cell-type_perc_dist.png"), p_perc, height = 10, width = 6, dpi = "print")
ggsave(here(output, "cell-type_perc_dist.pdf"), p_perc, height = 10, width = 6, dpi = "print")
