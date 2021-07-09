#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Plot genotpe PCA
# author: Jose Alquicira Hernandez
# date: 2021-03-27
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

#   local

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("uwot")

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-03-27", "population_umap")

#   ____________________________________________________________________________
#   Import data                                                             ####

ped <- read_tsv("data/genotype_pca/20130606_g1k.ped")
ev <- read_delim("data/genotype_pca/plink.eigenvec", 
                 delim = " ", 
                 col_names = c("Family ID", "Individual ID", "PC" %p% 1:20))


pop_codes <- read_tsv("Population Code	Population Description	Super Population Code	Sequence Data Available	Alignment Data Available	Variant Data Available
CHB	Han Chinese in Beijing, China	EAS	1	1	1
JPT	Japanese in Tokyo, Japan	EAS	1	1	1
CHS	Southern Han Chinese	EAS	1	1	1
CDX	Chinese Dai in Xishuangbanna, China	EAS	1	1	1
KHV	Kinh in Ho Chi Minh City, Vietnam	EAS	1	1	1
CEU	Utah Residents (CEPH) with Northern and Western European Ancestry	EUR	1	1	1
TSI	Toscani in Italia	EUR	1	1	1
FIN	Finnish in Finland	EUR	1	1	1
GBR	British in England and Scotland	EUR	1	1	1
IBS	Iberian Population in Spain	EUR	1	1	1
YRI	Yoruba in Ibadan, Nigeria	AFR	1	1	1
LWK	Luhya in Webuye, Kenya	AFR	1	1	1
GWD	Gambian in Western Divisions in the Gambia	AFR	1	1	1
MSL	Mende in Sierra Leone	AFR	1	1	1
ESN	Esan in Nigeria	AFR	1	1	1
ASW	Americans of African Ancestry in SW USA	AFR	1	1	1
ACB	African Caribbeans in Barbados	AFR	1	1	1
MXL	Mexican Ancestry from Los Angeles USA	AMR	1	1	1
PUR	Puerto Ricans from Puerto Rico	AMR	1	1	1
CLM	Colombians from Medellin, Colombia	AMR	1	1	1
PEL	Peruvians from Lima, Peru	AMR	1	1	1
GIH	Gujarati Indian from Houston, Texas	SAS	1	1	1
PJL	Punjabi from Lahore, Pakistan	SAS	1	1	1
BEB	Bengali from Bangladesh	SAS	1	1	1
STU	Sri Lankan Tamil from the UK	SAS	1	1	1
ITU	Indian Telugu from the UK	SAS	1	1	1
")

pop_codes <- pop_codes %>% select(Population = `Population Code`,  
                                  `Super Population Code`)



data <- left_join(ev, ped, by = "Individual ID")
data <- left_join(data, pop_codes, by = "Population")


data <- data %>% 
  mutate(Population = coalesce(Population, "OneK1K")) %>% 
  mutate(`Super Population Code` = coalesce(`Super Population Code`, "OneK1K"))


pop_order <- data$`Super Population Code` %>% unique() %>% sort()

pop_order[5] <- pop_order[6]
pop_order[6] <- "OneK1K"


data$`Super Population Code` <- factor(data$`Super Population Code`, pop_order)
data <- data %>% arrange(`Super Population Code`)



inds <- read.table("data/genotype_data/grm_subset.fam", 
                   sep = " ", 
                   col.names = c("famid", "indid", "paternal", 
                                 "maternal", "sex", "pheno"))

inds$id <- inds$famid %p% "_" %p% inds$indid


md <- readRDS("data/metadata/metadata.RDS")
md$individual[which(md$individual == "870_871" & md$latent == "b1")] <- "966_967"
inds <- md$individual %>% unique()

i <- (data$`Super Population Code` == "OneK1K") & (data$`Individual ID` %in% inds)
j <- data$`Super Population Code` != "OneK1K"

z <- i | j


data$status <- if_else(z, "final", "outlier")

data %>% 
  filter(`Super Population Code` == "OneK1K") %>% 
  nrow()

data %>% 
  filter(`Super Population Code` == "OneK1K", status == "final") %>% 
  nrow()

data <- data %>% filter(status == "final")

#   ____________________________________________________________________________
#   Plot data                                                               ####


p <- ggplot(data) +
  aes(x = PC1, 
      y = PC2, 
      color = `Super Population Code`, 
      shape = `Super Population Code`) +
  geom_point(alpha = 0.8, size = 3) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here(output, "population_pca.png"), p, width = 6.5, height = 5.5, dpi = "print")

# lims <- data %>% 
#   filter(`Super Population Code` == "OneK1K") %>% 
#   select(PC1, PC2) %>% 
#   summarise_all(list(mean = mean, sd = sd))
# p +
#   geom_vline(xintercept = lims$PC1_mean - 6*lims$PC1_sd) +
#   geom_vline(xintercept = lims$PC1_mean + 6*lims$PC1_sd) +
#   geom_hline(yintercept = lims$PC2_mean + 6*lims$PC2_sd) +
#   geom_hline(yintercept = lims$PC2_mean - 6*lims$PC2_sd)



#   ____________________________________________________________________________
#   Plot UMAP                                                               ####


pca <- data %>% 
  select(starts_with("PC")) %>% 
  as.data.frame() %>% 
  as.matrix()

umap_res  <- uwot::tumap(pca[,1:4], metric = "cosine", n_neighbors = 30, init = "spca")

umap_res %>% 
  as.data.frame() %>% 
  set_names("UMAP" %p% 1:2) %>% 
  mutate(pop = data$`Super Population Code`) -> umap_res

pal <- c("AFR" = "#63BA94", 
         "AMR" = "#EB8E66", 
         "EAS" = "#929FCB", 
         "EUR" = "#DA8CBC", 
         "SAS" = "#ACD161", 
         "OneK1K" = "#6FC6F5")


umap_res %>% 
  ggplot() +
  aes(x = UMAP1, 
      y = UMAP2, 
      color = pop) +
  geom_point(size = 1) +
  #geom_point(size = 1.5, shape = 21, color = "gray30", stroke = 0.01) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_color_manual(values = alpha(pal, 1)) +
  #scale_fill_manual(values = alpha(pal, 1)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Population",  
                              override.aes = list(size = 5)), 
         shape = guide_legend(title = "Population")) -> p

p

#   ____________________________________________________________________________
#   Export data                                                             ####


ggsave(here(output, "population_umap.png"), p, width = 6.5, height = 5.5, dpi = "print")
ggsave(here(output, "population_umap.pdf"), p, width = 6.5, height = 5.5, dpi = "print")


#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))

