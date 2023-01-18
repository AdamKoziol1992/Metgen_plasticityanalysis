setwd('~/Documents/GitHub/Metgen_plasticityanalysis/')
library(data.table)
library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(ggcorrplot)
library(factoextra)
library(ape)
library(logisticPCA)
library(vegan)
library(psych)
library(tidyverse)
library(janitor)
library(ggfortify)

`%!in%` = Negate(`%in%`)

# Load Apodemus module data
Apodemus_modules=fread("./distillR//outputsR/distilled_table_functions_AS.csv") %>% 
  column_to_rownames('V1')

# Load taxonomy
taxonomy <- fread("./data/taxa.csv", header = T) %>% 
  column_to_rownames('user_genome')

#Load Apodemus MAG quality data
AS_TPM <- read.csv('./data/AS_TPM.csv', header = T, row.names=1)

#Subset the good bins only
apodemus_damma <- Apodemus_modules[rownames(AS_TPM),] %>%
  mutate_all(as.numeric)

taxonomy_ordered <- taxonomy[rownames(apodemus_damma),] %>%
  as.data.frame() %>%
  rownames_to_column('bins') %>%
  mutate(Phylum = str_remove(Phylum, '_[A-B]')) %>% 
  select(-V1)

# PCA
module_pca=prcomp(apodemus_damma,scale = TRUE)
# variance explained by dimensions
fviz_screeplot(module_pca)

AS_PCA <- autoplot(module_pca, data = taxonomy_ordered, colour = 'Phylum', loadings = T, loadings.label = T) +
  theme_bw() +
  scale_colour_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#5ad45a", "#8be04e", "#ebdc78"))
AS_PCA
#Crocidura russula PCA

Crocidura_modules <- read.csv('./distillR/outputsR/distilled_table_functions_CR.csv',  header = T, row.names = 1)
CR_TPM <- read.csv('./data/TPM_CR.csv', row.names = 1)

#Filter only good MAGs
Crocidura_modules <- Crocidura_modules[rownames(CR_TPM),] %>%
  as.data.frame()

taxonomy_CR <- taxonomy[rownames(CR_TPM),] %>%
  as.data.frame() %>%
  mutate(Phylum = str_remove(Phylum, '_[A-B]'))
  
taxonomyCR_ordered <- taxonomy_CR[rownames(Crocidura_modules),]

# PCA
module_pca=prcomp(Crocidura_modules,scale = TRUE)
# variance explained by dimensions
fviz_screeplot(module_pca)

CR_PCA <- autoplot(module_pca, data = taxonomyCR_ordered, colour = 'Phylum', loadings = T, loadings.label = T) +
  theme_bw() +
  scale_colour_manual(values=c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#ffa300",  "#b3d4ff", "#8be04e", "#9b19f5"))

combined_PCA <- ggpubr::ggarrange(AS_PCA, CR_PCA, labels=c('A', 'B'), nrow = 2)
combined_PCA

ggsave('./PCA/PCA_combined.png', width = 10, height = 7)
ggsave('./PCA/PCA_combined.tiff', width = 10, height = 7, compression = 'lzw')
