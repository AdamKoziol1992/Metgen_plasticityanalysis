library(Hmsc)
library(colorspace)
library(corrplot)
library(gridExtra)
library(grid)
library(tidyverse)
setwd("~/Documents/GitHub/Metgen_plasticityanalysis/Hmsc/")
localDir = "~/Documents/GitHub/Metgen_plasticityanalysis/Hmsc/"
ModelDir = file.path(localDir, "models")
PanelDir = file.path(localDir, "panels")
DataDir = file.path(localDir, "data")

TH=2
thin = 10
samples = 250
nChains = 4

modelnames=c("CR") #"AS",

## Phylogenetic and functional structure on microbiome response to time
## ***********
load('./models/model_CR_TH2_thin_10_samples_250_chains_4.Rdata')
m_CR=m

i=2
load('./models/model_AS_TH2_thin_10_samples_250_chains_4.Rdata')
m_AS=m

## The effect of traits on microbiome response to time
VP_AS=computeVariancePartitioning(m_AS)
VP_AS$R2T
VP_CR=computeVariancePartitioning(m_CR)
VP_CR$R2T


## Phylogenetic structure on microbiome response to time (beyond the effect of traits)
mpost_AS=convertToCodaObject(m_AS)
quantile(unlist(mpost_AS$Rho),probs=c(.05,.5,.95))

mpost_CR=convertToCodaObject(m_CR)
quantile(unlist(mpost_CR$Rho),probs=c(.05,.5,.95))


beta_post_AS_beta=getPostEstimate(m_AS,parName = "Beta")
beta_post_CR_beta=getPostEstimate(m_CR,parName = "Beta")
plotBeta(m_AS,beta_post_AS_beta,plotTree = TRUE,supportLevel = 0.9)
plotBeta(m_CR,beta_post_CR_beta,plotTree = TRUE,supportLevel = 0.9)

## Functional dynamics

# AS
Gradient_Time_AS=constructGradient(m_AS,focalVariable = "Sampling.time",non.focalVariables = 1)
predY_Sampling.time_AS = predict(m_AS, Gradient = Gradient_Time_AS, expected = TRUE)

Polysaccharide.degradation_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                       showData = FALSE,main = "Polysaccharide degradation",index = 2,ylabel = "")
Sugar.degradation_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                          showData = FALSE,main = "Sugar degradation",index = 3,ylabel = "")
Lipid.degradation_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                               showData = FALSE,main = "Lipid degradation",index = 4,ylabel = "")
Protein.degradation_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                           showData = FALSE,main = "Protein degradation",index = 5,ylabel = "")
Mucin.degradation_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                showData = FALSE,main = "Mucin degradation",index = 6,ylabel = "")
SCFA.production_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                             showData = FALSE,main = "SCFA production",index = 7,ylabel = "")
Organic.anion.production_AS <- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                   showData = FALSE,main = "Organic anion production",index = 8,ylabel = "")
Secondary.bile.acid.production_AS<- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                     showData = FALSE,main = "Secondary.bile.acid.production",index = 9,ylabel = "")
Amino.acid.production_AS<- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                                     showData = FALSE,main = "Amino.acid.production",index = 10,ylabel = "")
Amino.acid.derivative.production_AS<- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                        showData = FALSE,main = "Amino.acid.derivative.production",index = 11,ylabel = "")
Vitamin.production_AS<- plotGradient(m_AS, Gradient_Time_AS, pred=predY_Sampling.time_AS, measure="T",q = c(0.05, 0.5, 0.95), 
                                        showData = FALSE,main = "Vitamin.production",index = 12,ylabel = "")


F1 <- as.data.frame(Polysaccharide.degradation_AS[[1]]) %>%
  mutate(Function = rep('Polysaccharide.degradation_AS'))
F2 <- as.data.frame(Sugar.degradation_AS[[1]]) %>%
  mutate(Function = rep('Sugar.degradation_AS'))
F3 <- as.data.frame(Lipid.degradation_AS[[1]]) %>%
  mutate(Function = rep('Lipid.degradation_AS'))
F4 <- as.data.frame(Protein.degradation_AS[[1]]) %>%
  mutate(Function = rep('Protein.degradation_AS'))
F5 <- as.data.frame(Mucin.degradation_AS[[1]]) %>%
  mutate(Function = rep('Mucin.degradation_AS'))
F6 <- as.data.frame(SCFA.production_AS[[1]]) %>%
  mutate(Function = rep('SCFA.production_AS'))
F7 <- as.data.frame(Organic.anion.production_AS[[1]]) %>%
  mutate(Function = rep('Organic.anion.production_AS'))
F8 <- as.data.frame(Secondary.bile.acid.production_AS[[1]]) %>%
  mutate(Function = rep('Secondary.bile.acid.production_AS'))
F9 <- as.data.frame(Amino.acid.production_AS[[1]]) %>%
  mutate(Function = rep('Amino.acid.production_AS'))
F10 <- as.data.frame(Amino.acid.derivative.production_AS[[1]]) %>%
  mutate(Function = rep('Amino.acid.derivative.production_AS'))
F11 <- as.data.frame(Vitamin.production_AS[[1]]) %>%
  mutate(Function = rep('Vitamin.production_AS'))

combined_functionabarplot_AS <- rbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11) %>%
  as.data.frame()
write.csv(combined_functionabarplot_AS, './functionalbarplots_AS.csv')
write.csv(combined_functionabarplot_AS, '../data/functionalbarplots_AS.csv')

windows(10)
df <- map(reduce(Polysaccharide.degradation_AS[[1]],
             Sugar.degradation_AS[[1]],
             Lipid.degradation_AS[[1]],
             Mucin.degradation_AS[[1]],
             SCFA.production_AS[[1]],
             Organic.anion.production_AS[[1]],
             Amino.acid.production_AS[[1]],
             Amino.acid.derivative.production_AS[[1]],
             Vitamin.production_AS[[1]]
))


Gradient_Time_CR=constructGradient(m_CR,focalVariable = "Sampling.time")
predY_Sampling.time_CR = predict(m_CR,Gradient = Gradient_Time_CR, expected = TRUE)

plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
             showData = FALSE)

F1 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                              showData = FALSE,main = "Amino acid biosynthesis",index = 2,ylabel = "")
F2 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Amino acid derivative biosynthesis",index = 3,ylabel = "")
F3 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "SCFA biosynthesis",index = 4,ylabel = "")
F4 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                       showData = FALSE,main = "Secondary bile acid biosynthesis",index = 5,ylabel = "")
F5 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Organic anion biosynthesis",index = 6,ylabel = "")
F6 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                   showData = FALSE,main = "Vitamin biosynthesis",index = 7,ylabel = "")
F7 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                            showData = FALSE,main = "Lipid degradation",index = 8,ylabel = "")
F8 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                                 showData = FALSE,main = "Polysaccharide degradation",index = 9,ylabel = "")
F9<- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                        showData = FALSE,main = "Protein degradation",index = 10,ylabel = "")
F10 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                                   showData = FALSE,main = "Sugar degradation",index = 11,ylabel = "")
F11 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Nitrogen compound degradation",index = 12,ylabel = "")
F12 <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                    showData = FALSE,main = "(Poly)amine degradation",index = 13,ylabel = "")

F1 <- as.data.frame(F1[[1]]) %>%
  mutate(Function = rep('Amino acid biosynthesis'))
F2 <- as.data.frame(F2[[1]]) %>%
  mutate(Function = rep('Amino acid derivative biosynthesis'))
F3 <- as.data.frame(F3[[1]]) %>%
  mutate(Function = rep('SCFA biosynthesis'))
F4 <- as.data.frame(F4[[1]]) %>%
  mutate(Function = rep('Secondary bile acid biosynthesis'))
F5 <- as.data.frame(F5[[1]]) %>%
  mutate(Function = rep('Organic anion biosynthesis'))
F6 <- as.data.frame(F6[[1]]) %>%
  mutate(Function = rep('Vitamin biosynthesis'))
F7 <- as.data.frame(F7[[1]]) %>%
  mutate(Function = rep('Lipid degradation'))
F8 <- as.data.frame(F8[[1]]) %>%
  mutate(Function = rep('Polysaccharide degradation'))
F9 <- as.data.frame(F9[[1]]) %>%
  mutate(Function = rep('Protein degradation'))
F10 <- as.data.frame(F10[[1]]) %>%
  mutate(Function = rep('Sugar degradation'))
F11 <- as.data.frame(F11[[1]]) %>%
  mutate(Function = rep('Nitrogen compound degradation'))
F12 <- as.data.frame(F12[[1]]) %>%
  mutate(Function = rep('(Poly)amine degradation'))

combined_functionabarplot_CR <- rbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12) %>%
  as.data.frame() %>% 
  mutate(Function = str_replace_all(Function, ' ', '_'),
         Function = str_c(Function, '_CR'))
write.csv(combined_functionabarplot_CR, '../data/functional_barplot_CR')


## MAG-level responses
# CR

Gradient_Time_CR=constructGradient(m_CR,focalVariable = "Sampling.time",non.focalVariables = 1)
predY_Sampling.time_CR = predict(m_CR,Gradient = Gradient_Time_CR, expected = TRUE)
library(abind)
predY_CR_time=abind(predY_Sampling.time_CR,along=3)
predY_CR_time_median=apply(predY_CR_time, 1:2, median)

# Subset the mags that increased in relative abundance in acclimation period and 
# then again in heat treatment
index=which(predY_CR_time_median[1,]<predY_CR_time_median[2,]&predY_CR_time_median[2,]<predY_CR_time_median[3,])

# Manually exclude MAGs decreasing during heat after the initial increase in acclimation
# (based on the visualization of below plots)

# List of MAGs increasing from baseline to acclimation and from acclimation to heat.
names(index)

plot_list=list()
for(i in 1:length(index)){
  temp_index=index[i]
  p=plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="Y",q = c(0.05, 0.5, 0.95), 
               showData = FALSE,index=temp_index)
  plot_list[[i]]=p
}
c1=plot_list[[1]];c2=plot_list[[2]];c3=plot_list[[3]];c4=plot_list[[4]];c5=plot_list[[5]];
c6=plot_list[[6]];c7=plot_list[[7]];c8=plot_list[[8]];c9=plot_list[[9]];c10=plot_list[[10]];
c11=plot_list[[11]];c12=plot_list[[12]];c13=plot_list[[13]];c14=plot_list[[14]];
c15=plot_list[[15]];c16=plot_list[[16]];c17=plot_list[[17]]

library(gridExtra)
library(grid)
grid.arrange(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,
                           ncol=4,nrow=5)

# Subset the mags that decreased in relative abundance in acclimation period and 
# then again in heat treatment
index=which(predY_CR_time_median[1,]>predY_CR_time_median[2,]&predY_CR_time_median[2,]>predY_CR_time_median[3,])

# List of MAGs increasing from baseline to acclimation and from acclimation to heat.
names(index)

plot_list=list()
for(i in 1:length(index)){
  temp_index=index[i]
  p=plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="Y",q = c(0.05, 0.5, 0.95), 
                 showData = FALSE,index=temp_index)
  plot_list[[i]]=p
}
c1=plot_list[[1]];c2=plot_list[[2]];c3=plot_list[[3]];c4=plot_list[[4]];c5=plot_list[[5]];
c6=plot_list[[6]];c7=plot_list[[7]];c8=plot_list[[8]];c9=plot_list[[9]];c10=plot_list[[10]];
c11=plot_list[[11]];c12=plot_list[[12]];c13=plot_list[[13]];c14=plot_list[[14]];
c15=plot_list[[15]];c16=plot_list[[16]]

library(gridExtra)
library(grid)
grid.arrange(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,
             ncol=4,nrow=4)

plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="S",q = c(0.05, 0.5, 0.95), 
             showData = FALSE)

Polysaccharide.degradation_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                              showData = FALSE,main = "Polysaccharide degradation",index = 2,ylabel = "")
Sugar.degradation_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Sugar degradation",index = 3,ylabel = "")
Lipid.degradation_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Lipid degradation",index = 4,ylabel = "")
Mucin.degradation_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Mucin degradation",index = 5,ylabel = "")
SCFA.production_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                   showData = FALSE,main = "SCFA production",index = 6,ylabel = "")
Organic.anion.production_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                            showData = FALSE,main = "Organic anion production",index = 7,ylabel = "")
Amino.acid.production_CR <- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                         showData = FALSE,main = "Amino acid production",index = 8,ylabel = "")
Amino.acid.derivative.production_CR<- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                                   showData = FALSE,main = "Amino acid derivative production",index = 9,ylabel = "")
Vitamin.production_CR<- plotGradient(m_CR, Gradient_Time_CR, pred=predY_Sampling.time_CR, measure="T",q = c(0.05, 0.5, 0.95), 
                                     showData = FALSE,main = "Vitamin production",index = 10,ylabel = "")
grid.arrange(Polysaccharide.degradation_CR,
             Sugar.degradation_CR,
             Lipid.degradation_CR,
             Mucin.degradation_CR,
             SCFA.production_CR,
             Organic.anion.production_CR,
             Amino.acid.production_CR,
             Amino.acid.derivative.production_CR,
             Vitamin.production_CR,
             ncol=3,nrow=3)



