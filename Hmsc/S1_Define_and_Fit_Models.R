library(Hmsc)
library(tidyverse)
setwd <- "~/Documents/GitHub/Metgen_plasticityanalysis/Hmsc/"

# Y matrix
TH=2 #TH=2 for bigger data; TH=1 for biggest data
if(TH==3){
  thresh_abu=0.05
}
if(TH==2){
  thresh_abu=0.01
}
if(TH==1){
  thresh_abu=0.005
}

## Load and prepare AS data
load(file.path(getwd(),"data/allData_AS.R"))
YData_AS=YData
XData_AS=XData
SData_AS=SData
PData_AS=PData
TrData_AS=TrData

#### AS Ydata based on treatments
YData_AS=YData_AS[,apply(YData_AS,2,max)>thresh_abu]
dim(YData_AS)
YData_AS=log1p(YData_AS)

#SData
SData_AS <- SData_AS %>%
  mutate_all(as.factor)
###Xdata
XData_AS <- XData_AS %>%
  as.data.frame
rownames(SData_AS)=rownames(XData_AS)

## Load and prepare CR data
load(file.path(getwd(),"data/allData_CR.R"))
YData_CR=YData
XData_CR=XData
SData_CR=SData
PData_CR=PData
TrData_CR=TrData
## Define matrices of the Hmsc model
## *********************************

#### CR Ydata based on treatments
YData_CR=YData_CR[,apply(YData_CR,2,max)>thresh_abu]
dim(YData_CR)
YData_CR=log1p(YData_CR)
#SData
SData_CR <- SData_CR %>% mutate_all(as.factor)

###Xdata
XData_CR <- XData_CR %>%
  as.data.frame

#rownames(SData_CR)=rownames(XData_CR)


## Define formulas of the Hmsc model
## *********************************

# X
XFormula=~Sampling.time

# Tr
TrFormula=~B03+
  B04+
  B01+
  B05+
  B02+
  B06+
  D03+
  D01+
  D04+
  D02+
  D08+
  D05
  
# StudyDesign
rL.Pen_AS = HmscRandomLevel(units = levels(SData_AS[,1]))
rL.Individual_AS = HmscRandomLevel(units = levels(SData_AS[,2]))
rL.Pen_CR = HmscRandomLevel(units = levels(SData_CR[,1]))
rL.Individual_CR = HmscRandomLevel(units = levels(SData_CR[,2]))

## Define Hmsc models ###### I included the phylogenetic tree but don't have it in this model
## ******************
m_AS = Hmsc(Y=YData_AS,XData = XData_AS, XFormula = XFormula, studyDesign = SData_AS,phyloTree = PData_AS, 
         ranLevels = list("cage"=rL.Pen_AS,"individual_id"=rL.Individual_AS),
         TrData = TrData_AS, TrFormula = TrFormula, distr = "normal",YScale = TRUE)

m_CR = Hmsc(Y=YData_CR,XData = XData_CR, XFormula = XFormula, studyDesign = SData_CR,phyloTree = PData_CR, 
            ranLevels = list("cage"=rL.Pen_CR,"individual_id"=rL.Individual_CR),
            TrData = TrData_CR, TrFormula = TrFormula, distr = "normal",YScale = TRUE)

#models=list(m_AS,m_CR)
#modelnames=c("AS","CR")

models=list(m_CR)
modelnames=c("CR")

## Fit models
## **********

samples_list = c(10, 250)
thin_list = c(1, 10)
nChains = 4

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  for(i in 1:length(models)){
    m = sampleMcmc(models[[i]], samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),models[[i]]$nr),
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains,nParallel = nChains)
    filename=paste("model_",modelnames[i],"_TH",as.character(TH),"_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
    save(m,file=file.path(getwd(),'models',filename))
  }
}

models=list(m_AS)
modelnames=c("AS")

## Fit models
## **********

samples_list = c(10, 250)
thin_list = c(1, 10)
nChains = 4

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  for(i in 1:length(models)){
    m = sampleMcmc(models[[i]], samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),models[[i]]$nr),
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains,nParallel = nChains)
    filename=paste("model_",modelnames[i],"_TH",as.character(TH),"_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
    save(m,file=file.path(getwd(),'models',filename))
  }
}
