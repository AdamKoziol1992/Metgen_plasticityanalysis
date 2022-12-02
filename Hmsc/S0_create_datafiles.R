## 1. Load and prepare datasets
## ****************************
setwd('~/Documents/GitHub/Metgen_plasticityanalysis/Hmsc')
library(tidyverse)
library(janitor)
library(distillR)
library(ape)
library(phytools)

f_table <- pathway_db %>%
  mutate(Function = str_replace_all(Function, ' ', '_')) ##Indolic compounds has been merged with Aromatic_compounds

`%!in%` = Negate(`%in%`)

Functions_to_keep <- f_table %>%
  select(Function) %>%
  filter(Function %!in% c('Protein_degradation', 'Antibiotic_biosynthesis', 'Metallophore_biosynthesis', 'Nucleic_acid_biosynthesis')) %>%
  left_join(f_table[,c('Function', 'Code_function')], by = 'Function') %>% 
  mutate(Function = str_replace(Function, 'Indolic_compound_biosynthesis', 'Aromatic_compound_biosynthesis'),
         Code_function = str_replace(Code_function, 'B05', 'B08')) %>% 
  unique %>% 
  arrange(desc(Code_function))


### 1.1. Load MAG count data ###
MAG_counts_AS <- read.csv("../data/AS_TPM.csv",sep = ",", row.names = 1) %>%
  select(!contains('FAST')) %>%
  as.data.frame

MAG_counts_CR <- read.csv('../data/TPM_CR.csv', sep = ',', row.names = 1)

final_mag_counts_as <- MAG_counts_AS
final_mag_counts_cr <- MAG_counts_CR

### 1.2. Load metadata and reorder the rows ###
metadata_as=read.csv('../data/metadata.csv', header = T, row.names = 1) %>%
  filter(species == 'AS') %>%
  column_to_rownames('sample_id') %>% 
  .[colnames(final_mag_counts_as),]

metadata_cr=read.csv('../data/metadata.csv', header = T, row.names = 1) %>%
  filter(species == 'CR') %>%
  column_to_rownames('sample_id') %>% 
  .[colnames(final_mag_counts_cr),]


# Edit the metadata to include time as a factor
metadata_as <- metadata_as %>% 
  mutate(treatment = factor(treatment, levels = c('Day1', 'Acc', 'Heat', 'Cold', 'Diet')),
         time_point = case_when(grepl('^Day1', metadata_as$treatment) ~ 0,
                                grepl('^Acc', metadata_as$treatment) ~ 7,
                                grepl('^Heat', metadata_as$treatment) ~ 19,
                                grepl('^Cold', metadata_as$treatment) ~ 33,
                                grepl('^Diet', metadata_as$treatment) ~ 45)) %>% 
  mutate_at(c('time_point', 'individual_id', 'cage', 'sex'), as.factor)

metadata_cr <- metadata_cr %>% 
  mutate(treatment = factor(treatment, levels = c('Day1', 'Acc', 'Heat', 'Cold', 'Diet')),
         time_point = case_when(grepl('^Day1', metadata_cr$treatment) ~ 0,
                                grepl('^Acc', metadata_cr$treatment) ~ 7,
                                grepl('^Heat', metadata_cr$treatment) ~ 19,
                                grepl('^Cold', metadata_cr$treatment) ~ 33,
                                grepl('^Diet', metadata_cr$treatment) ~ 45)) %>% 
  mutate_at(c('time_point', 'individual_id', 'cage', 'sex'), as.factor)

# 1.4 Load phylogenetic trees
p_tree_AS <- read.tree(file='../data/gtdbtk.bac120.classify_AS.tree') %>%
  keep.tip(rownames(final_mag_counts_as)) %>%
  force.ultrametric(p_tree_AS, method="extend")

p_tree_CR <- read.tree(file='../data/gtdbtk.bac120.classify_CR.tree') %>%
  keep.tip(rownames(final_mag_counts_cr)) %>%
  force.ultrametric(p_tree_AS, method="extend")

##Re order the species in count matrix to that of the phylogenetic tree
final_mag_counts_as <- final_mag_counts_as %>%
  .[p_tree_AS$tip.label,]

final_mag_counts_cr <- final_mag_counts_cr %>%
  .[p_tree_CR$tip.label,]

### 1.3. Load functional annotation table ###

MAG_annot_as=read.csv("../distillR/outputsR/distilled_table_functions_AS.csv",sep = ",", row.names = 1) %>% 
  .[rownames(final_mag_counts_as),] %>%
  .[,Functions_to_keep$Code_function] %>% 
  mutate_all(as.numeric)
MAG_annot_cr=read.csv("../distillR/outputsR/distilled_table_functions_CR.csv",sep = ",", row.names = 1) %>% 
  .[rownames(final_mag_counts_cr),] %>%
  .[,Functions_to_keep$Code_function] %>% 
  mutate_all(as.numeric)

paste("MAGs in MAG_weighted_ca and MAG_annot tables in same order:",
      mean(rownames(MAG_annot_as)==rownames(final_mag_counts_as))==1
)

paste("MAGs in MAG_weighted_ca and MAG_annot tables in same order:",
      mean(rownames(MAG_annot_cr)==rownames(final_mag_counts_cr))==1
)

## Save data matrices to be used in the Hmsc modeling
## **************************************************

# Y matrix
YData=t(final_mag_counts_as)

# X matrix
XData=data.frame(Sampling.time=metadata_as$time_point, mouse=metadata_as$individual_id)

# S matrix
SData=data.frame(cage=metadata_as$cage,individual_id=metadata_as$individual_id,sample_type=metadata_as$treatment)

# Tr matrix
TrData=data.frame(MAG_annot_as)
rownames(TrData)=rownames(final_mag_counts_as)
rownames(TrData)==colnames(YData)

# P matrix
PData <- p_tree_AS

# Save the datafiles
save(YData,XData,SData,TrData,PData,
     file="../data/allData_AS.R")

## Save data matrices to be used in the Hmsc modeling
## Crocidura russula
## **************************************************

# Y matrix
YData=t(final_mag_counts_cr)

# X matrix
XData=data.frame(Sampling.time=metadata_cr$time_point, shrew=metadata_cr$individual_id)

# S matrix
SData=data.frame(cage=metadata_cr$cage,individual_id=metadata_cr$individual_id,sample_type=metadata_cr$treatment)

# Tr matrix
TrData=data.frame(MAG_annot_cr)
rownames(TrData)=rownames(final_mag_counts_cr)
rownames(TrData)==colnames(YData)

# P matrix
PData <- p_tree_CR

#Save the data files
save(YData,XData,SData,TrData,PData,
     file="../data/allData_CR.R")

