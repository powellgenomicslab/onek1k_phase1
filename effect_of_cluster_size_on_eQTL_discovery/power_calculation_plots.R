##############################################################################
# Script information                                                      
# Title: Investigate results of cluster size on eQTL discovery (Fig. S17)
# Author: Seyhan Yazar
# Date: 2021-04-08
# Description: None
##############################################################################


# Upload libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(magrittr)

wd <- '/onek1k/effect_of_cluster_size_on_eQTL_discovery/percent_specific_results'
setwd(wd)
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.tsv")
dataset <- lapply(files, function(x) { fread(x)})
# Get file names and add to the list
filenames <- sub("*.tsv", "", files)
names(dataset) <- filenames
df_eqtl <- bind_rows(dataset, .id="df.name")
df_eqtl$percent <- gsub("_.*","", df_eqtl$df.name)
df_eqtl$eSNP_rank <- gsub(".*_","",df_eqtl$df.name)

df_100 <- fread("/onek1k/single_cell_cis_eQTL_mapping/supplementary/OneK1K_eQTLs_Results_03022021.tsv")
df_100 <- df_100 %>% filter(cell_type=="CD4 NC")

pcount <- data.frame(df_eqtl %>% group_by(percent) %>% count())
pcount[nrow(pcount) + 1,] = list("percent100", 6473)
pcount$percent2 <- as.numeric(gsub('percent*',"", pcount$percent))
pcount <- pcount %>% arrange(-percent2)
pcount$total <- as.numeric(c("463528", "347651", "231763", "115878", "46352", "23187", "4648"))
pcount$n <- as.numeric(pcount$n)
pcount <- pcount[,-1]
pcount$celltype <- "CD4 NC"
pcount$ave <- round(pcount$total/982,0)

pointdata <- data.frame(table(df_100$cell_type))
colnames(pointdata) <- c("celltype", "n")
pointdata$total <- as.numeric(c("82068","48023","61786","463528","4065","205077","133482","34528","8690","38233","15166", "159820","9677","3625"))
pointdata$percent2 <- round((pointdata$total/463528)*100, 0)
pointdata <- pointdata[-4,]
pointdata <- pointdata %>% arrange(-percent2)
pointdata$n <- as.numeric(pointdata$n)

bothdf <- rbind(pcount, pointdata)

# Set order of plotting
corder <- c("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET","CD8 NC", "CD8 S100B", "NK", "NK R", "Plasma","B Mem", "B IN","Mono C", "Mono NC","DC")
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")

tol14rainbow2 <- data.frame(cbind(corder, tol14rainbow))  
colnames(tol14rainbow2) <- c("celltype","rainbow") 
tol14rainbow2$celltype <- as.character(tol14rainbow2$celltype)
tol14rainbow2$rainbow <- as.character(tol14rainbow2$rainbow)

pointdata$celltype <- as.character(pointdata$celltype)   
pointdata <- left_join(pointdata, tol14rainbow2, by="celltype")
     
###### New plot with percentages ######
colnames(pcount)[2] <- "percent"
pcount <- left_join(pcount,gene_no,by="percent")
pcount$percentage <- round((pcount$n/pcount$after_filtering)*100, 2)

setwd("/onek1k/single_cell_cis_eQTL_mapping/supplementary/gene_tested")
# Get all the result files (significant eQTls)
files <- list.files(path=".", pattern="*.txt")
dataset2 <- lapply(files, function(x) { fread(x)})
filenames <- sub("*.txt", "", files)
names(dataset2) <- filenames
df <- bind_rows(dataset2, .id="dfname")
genes_or <- df %>% group_by(dfname) %>% count(n())
genes_or$celltype <- c("B IN", "B Mem", "CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", "CD8 ET", "CD8 S100B", "DC", "Mono C", "Mono NC", "NK R", "NK", "Plasma")
genes_or <- genes_or %>% select (celltype, n)

pointdata <- left_join(pointdata, genes_or, by="celltype")
pointdata$percentage <- round((pointdata$n.x/pointdata$n.y)*100, 2)

png("/onek1k/effect_of_cluster_size_on_eQTL_discovery/power_plot_210408.png", width=12, height=8)
ggplot() +
    theme_classic(14) +
    scale_x_continuous(lim=c(0,100), breaks=c(1,5,10,25,50,75,100),
    labels=c('1%\n(n=5)','5%\n(n=24)', '10%\n(n=47)', 
    '25%\n(n=118)','50%\n(n=236)', '75%\n(n=354)', '100%\n(n=472)')) +
    geom_point(data=pcount, aes(x=percent, y=percentage), col="#882E72", size=3) +
    geom_line(data=pcount, aes(x=percent, y=percentage), col="#882E72", size=1) +
    geom_point(data=pointdata, aes(x=percent2, y=percentage), col=pointdata$rainbow, size=3) +
    geom_text_repel(data=pointdata, aes(x=percent2, y=percentage, label=celltype), force=3) +
    labs(y="Percentage of independent eQTLs over genes tested in CD4 NC cells", x="")
dev.off()



