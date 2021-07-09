#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Identify OneK1K dynamic eQTLs in B cells
# author: Jose Alquicira Hernandez
# date: 2021-04-16
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# screen -S dyn-eqtl
# qrsh -N dyn-eqtl -l mem_requested=100G -pe smp 1 -q short.q
# conda activate r-4.0.3

#   ____________________________________________________________________________
#   Import libraries                                                        ####

# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("data.table")
library("furrr")
library("broom")
library("ComplexHeatmap")
library("ggvenn")


#   ____________________________________________________________________________
#   Future settings                                                         ####

options(future.globals.maxSize = 1024^3*100)
plan(multicore(workers = 2))

#   ____________________________________________________________________________
#   Set output                                                              ####

output <- set_output("2021-04-27", "dyn_eqtls")

#   ____________________________________________________________________________
#   Import gene expression data                                             ####


path <- file.path("/",
                  "directflow", 
                  "SCCGGroupShare", 
                  "projects", 
                  "SeyhanYazar",
                  "onek1k",
                  "science_revision_Dec20",
                  "pseudotime_analysis", 
                  "bcell_six_quantile")

exp_quantile_dirs <- list.files(path, full.names = TRUE, pattern = "Q")


quantile_names <- exp_quantile_dirs %>%  str_split("/") %>% map_chr(~.[length(.)])

exp_data <- future_map(exp_quantile_dirs, ~{
  
  message("Importing data from... \n", .x)
  
  files <- list.files(.x, pattern = "*log_residuals.tsv", full.names = TRUE)[1]
  mats <- fread(file = files, data.table = FALSE)
  mats <- mats %>% column_to_rownames("sampleid")
  mats
  
})


exp <- map(exp_data, colMeans)
exp <- map(exp, ~{ tibble(gene = names(.), residual = .)})
names(exp) <- quantile_names

exp <- bind_rows(exp, .id = "quantile")

#   ____________________________________________________________________________
#   Import data                                                             ####

path <- file.path("/",
                  "directflow", 
                  "SCCGGroupShare", 
                  "projects", 
                  "SeyhanYazar",
                  "onek1k",
                  "science_revision_Dec20",
                  "pseudotime_analysis",
                  "bcell_six_quantile",
                  "matrixeqtl_outputs")


quantile_dirs <- list.files(path, full.names = TRUE)

quantile_data <- future_map(quantile_dirs, ~{
  message("Importing data for", .x)
  files <- list.files(.x, pattern = "*.tsv", full.names = TRUE)
  chr <- files %>% str_split("/") %>% map(~.[length(.)]) %>% str_split("_") %>% 
    map(2)
  names(files) <- chr
  mats <- map(files, fread, 
              colClasses = c("character", "character", rep("numeric", 4)))
  rbindlist(mats, idcol = "chr")
})


names(quantile_data) <- quantile_dirs %>% str_split("/") %>% map_chr(~.[length(.)])


# Combine data across quantiles
data <- rbindlist(quantile_data, idcol = "quantile")


# Index data by snp-gene pair
data[,id :=  paste(SNP, gene, sep = "-")]

inicio("Indexing...")
setkey(data, id)
fin()
# ★  Elapsed time: 5.713 mins

inicio("Save eQTL aggregated data")
saveRDS(data, here(output, "eQTL.RDS"))
#data <- readRDS(here(output, "eQTL.RDS"))
fin()

n_eqtl <- data[FDR < 0.05, .N, by = "quantile"][order(quantile),]

# quantile     N
# 1:       Q1 29983
# 2:       Q2 25459
# 3:       Q3 26046
# 4:       Q4 29345
# 5:       Q5 29330
# 6:       Q6 42670

summary(n_eqtl[["N"]])


#   ____________________________________________________________________________
#   Import conditional analysis results                                     ####

path <- file.path("/",
                  "directflow", 
                  "SCCGGroupShare", 
                  "projects", 
                  "SeyhanYazar",
                  "onek1k",
                  "science_revision_Dec20",
                  "pseudotime_analysis",
                  "bcell_six_quantile",
                  "spearmans_output",
                  "pseudotime_conditional_analysis_results_210425.tsv")

ceqtl <- fread(file = path)
ceqtl <- ceqtl %>% mutate(id = paste(snpid, geneid, sep = "-"))


ceqtl[, .N, by = "eSNP_rank"][["N"]] %>% summary()


# Select SNP-gene pair info
ceqtl <- ceqtl %>% 
  select(id, snpid, geneid, rsid = ID) %>% 
  distinct()


# Subset conditionl eQTLs from MatrixeQTL results
data_sig <- data[ceqtl$id] %>% as_tibble()

# Add rsids
data_sig <- full_join(data_sig, ceqtl[, c("id", "rsid")], by = "id")

# Create explanatory numeric variable
data_sig$rank <- data_sig$quantile %>% str_remove("Q") %>% as.integer()

# Get SNP-gene pairs present in all quantiles
data_filter <- data_sig %>% 
  group_by(id) %>% 
  tally() %>% 
  filter(n == 6)

# Subset ubiquitous pairs 
data_sig <- data_sig[data_sig$id %in% data_filter$id, ]

# Nest data by SNP-gene pairs
data_sig <- data_sig %>% 
  group_by(rsid, gene, id) %>% 
  nest()


#   ____________________________________________________________________________
#   Identify dynamic eQTL                                                   ####

##  ............................................................................
##  Linear regression                                                       ####

data_sig <- data_sig %>% 
  mutate(linear_fit = map(data, ~{lm(beta ~ rank, .)}), 
         tidy_linear_fit = map(linear_fit, tidy), 
         glance_linear_fit = map(linear_fit, glance)
  )


lm_res <- data_sig %>% 
  unnest("glance_linear_fit") %>% 
  select(id, rsid, gene, r.squared, statistic, p.value) %>% 
  filter(p.value < 0.05) %>% 
  arrange(abs(statistic), desc(r.squared))

##  ............................................................................
##  Quadratic regression                                                    ####


data_sig <- data_sig %>% 
  mutate(quad_fit = map(data, ~{lm(beta ~ poly(rank, 2), .)}), 
         tidy_quad_fit = map(quad_fit, tidy),
         glance_quad_fit = map(quad_fit, glance)
  )

quad_res <- data_sig %>% 
  unnest("glance_quad_fit") %>% 
  select(id, rsid, gene, r.squared, statistic, p.value) %>% 
  filter(p.value < 0.05) %>% 
  arrange(abs(statistic), desc(r.squared))

#   ____________________________________________________________________________
#   Aggregate regression results                                            ####


data_sig$pair <- paste(data_sig$rsid, data_sig$gene, sep = "-")

res <- data_sig %>% 
  ungroup() %>% 
  select(id, pair, rsid, gene, data, 
         linear = glance_linear_fit, quadratic = glance_quad_fit) %>% 
  unnest("linear", names_sep = "_") %>% 
  unnest("quadratic", names_sep = "_") %>% 
  select(id, pair, rsid, gene,
         linear_statistic, linear_p.value, linear_r.squared,
         quadratic_statistic, quadratic_p.value, quadratic_r.squared) %>% 
  mutate(linear_significant = if_else(linear_p.value < 0.05, TRUE, FALSE),
         quadratic_significant = if_else(quadratic_p.value < 0.05, TRUE, FALSE),
         trajectory = if_else(linear_significant | quadratic_significant, TRUE, FALSE))




# Extract betas and transform to matrix format
betas <- data_sig %>% 
  ungroup() %>% 
  select(pair, data) %>% 
  unnest("data") %>% 
  pivot_wider(id_cols = "pair", 
              names_from = quantile,
              values_from = beta)


res <- full_join(res, betas, by = "pair")

write_tsv(res %>% select(-id), here(output, "linear_quad_eQTL.tsv"))

write_tsv(res %>% select(-id) %>% filter(trajectory), 
          here(output, "trajectory_eQTL.tsv"))


#   ____________________________________________________________________________
#   Plot intersection of detected trajectory SNP-gene pairs                 ####


lm_res_gene <- res %>% filter(linear_p.value < 0.05) %>% pull(gene) %>% unique()
quad_res_gene <- res %>% filter(quadratic_p.value < 0.05) %>% pull(gene) %>% unique()

res_genes <- list(lm = lm_res_gene, 
                  quad = quad_res_gene)


p <- ggvenn(
  res_genes, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(here(output, "genes_shared.png"), p + ggtitle("Genes"))


lm_assoc <- res %>% filter(linear_p.value < 0.05) %>% pull(id)
quad_assoc <- res %>% filter(quadratic_p.value < 0.05) %>% pull(id)

res_assoc <- list(lm = lm_assoc, 
                  quad = quad_assoc)

p <- ggvenn(
  res_assoc, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(here(output, "assoc_shared.png"), p + ggtitle("SNP-gene pairs"))



#   ____________________________________________________________________________
#   Plot beta estimates for significant results                             ####


linear_scaled_betas <- res %>% 
  ungroup() %>% 
  filter(linear_significant) %>% 
  select(pair, starts_with("Q", ignore.case = FALSE)) %>% 
  column_to_rownames("pair") %>% 
  apply(1, scale) %>% 
  t() %>% 
  `colnames<-`("Q" %p% 1:6) %>% 
  as.data.frame() %>% 
  rownames_to_column("pair") %>% 
  mutate(method = "Linear")


quadratic_scaled_betas <- res %>% 
  ungroup() %>% 
  filter(quadratic_significant) %>% 
  select(pair, starts_with("Q", ignore.case = FALSE)) %>% 
  column_to_rownames("pair") %>% 
  apply(1, scale) %>% 
  t() %>% 
  `colnames<-`("Q" %p% 1:6) %>% 
  as.data.frame() %>% 
  rownames_to_column("pair") %>% 
  mutate(method = "Quadratic")


shared_genes <- intersect(linear_scaled_betas$pair, quadratic_scaled_betas$pair)


all_scaled_betas <- bind_rows(linear_scaled_betas, 
            quadratic_scaled_betas %>% filter(!pair %in% shared_genes))


all_scaled_betas_mat <- all_scaled_betas %>% 
  select(pair, starts_with("Q", ignore.case = FALSE)) %>% 
  column_to_rownames("pair") %>% 
  apply(1, scale) %>% 
  t() %>% 
  `colnames<-`("Q" %p% 1:6)


all_scaled_betas$method <- factor(all_scaled_betas$method, 
                                  levels = rev(c("Linear", "Quadratic")))

# pdf(here(output, "all_scaled_betas_RdYlBu.pdf"), height = 9)
# Heatmap(all_scaled_betas_mat, 
#         col = RColorBrewer::brewer.pal(11, "RdYlBu"),
#         cluster_columns = FALSE, 
#         show_row_names = FALSE, 
#         row_split = all_scaled_betas$method, 
#         row_gap = unit(8, "mm"),
#         heatmap_legend_param = list(title = expression("Scaled "*beta*" values")))
# dev.off()





gwas <- fread(here("data", "b_cells_gwas", "B_cell_GWAS_dynamic_GWAS_7AI.csv"))
gwas[,pair := paste(ID, gene, sep = "-")]
setkey(gwas, pair)


gwas <- gwas[, .(pair, GWAS_Trait)]
gwas <- unique(gwas, by = c("pair", "GWAS_Trait"))
gwas <- gwas[,.(pair, GWAS_Trait)][, .N, by = "pair"]

lintersect(all_scaled_betas$pair,  gwas$pair)
setdiff(gwas$pair, all_scaled_betas$pair)

all_scaled_betas <- left_join(all_scaled_betas, gwas, by = "pair")



all_scaled_betas <- all_scaled_betas %>% mutate(gwas = if_else(is.na(N), "False", "True")) 


model <-  factor(all_scaled_betas$method, 
                 levels = c("Linear", "Quadratic"))

gwas <-  factor(all_scaled_betas$gwas, 
                 levels = c("True", "False"))
table(gwas)

ha <-  HeatmapAnnotation(Model = model,
                       col = list(Model = c(Linear = "#38618c", 
                                            Quadratic = "#ff5964")
                                  ),
                       annotation_label = ""
                       )



hb <-  HeatmapAnnotation(`Autoimmune GWAS` = gwas,
                         col = list(
                           `Autoimmune GWAS` = c(False = "#296A6C", 
                                             True = "#F89423")
                         )
)



pdf(here(output, "all_scaled_betas_RdYlBu_transp.pdf"), height = 3, width = 12)
Heatmap(all_scaled_betas_mat %>% t(),
        col = RColorBrewer::brewer.pal(11, "RdYlBu"),
        top_annotation = ha,
        bottom_annotation = hb,
        row_title = "Pseudotime",
        column_title = "eQTL",
        column_title_side = "bottom",
        cluster_rows  = FALSE, 
        #show_column_names = FALSE, 
        column_names_gp = grid::gpar(fontsize = 2),
        column_split = all_scaled_betas$method, 
        column_gap = unit(0, "mm"), 
        row_names_side = "left",
        heatmap_legend_param = list(title = expression("Scaled "*beta*" values")))
dev.off()


#   ____________________________________________________________________________
#   Plot trajectory eQTLs                                                   ####

res_exp <- res %>% 
  ungroup() %>% 
  filter(trajectory) %>% 
  select(id, pair, rsid, gene , starts_with("Q", ignore.case = FALSE)) %>% 
  pivot_longer("Q" %p% 1:6, names_to = "quantile", values_to = "beta") %>% 
  left_join(exp, by = c("quantile", "gene"))

res_exp$rank <- str_remove(res_exp$quantile, "Q") %>% as.integer()

res_exp <- split(res_exp, res_exp$pair) %>% map(as.data.frame)

dir.create(here(output, "plots"))
dir.create(here(output, "plots", "linear"))
dir.create(here(output, "plots", "quad"))


plot_pair <- function(x, name){
  
  p_quad <- ggplot(x, aes(rank, beta)) +
    geom_point(size = 3, color = "gray20") +
    stat_smooth(method = "lm", se = FALSE, color = "red", formula = y ~ poly(x, 2)) +
    xlab("Quantile") +
    ylab(expression(beta)) +
    ggtitle(name) +
    scale_x_continuous(breaks = 1:6) +
    theme_pub()
  ggsave(here(output, "plots",  "quad", name %p% ".pdf"), p_quad, height = 4.5, width = 6)
  
  
  p_lm <- ggplot(x, aes(rank, beta)) +
    geom_point(size = 3, color = "gray20") +
    stat_smooth(method = "lm", se = FALSE, color = "red") +
    xlab("Quantile") +
    ylab(expression(beta)) +
    ggtitle(name) +
    scale_x_continuous(breaks = 1:6) +
    theme_pub()
  ggsave(here(output, "plots",  "linear", name %p% ".pdf"), p_lm, height = 4.5, width = 6)
  
}


iwalk(res_exp, plot_pair)


#   ____________________________________________________________________________
#   Session info                                                            ####

print_session(here(output))


