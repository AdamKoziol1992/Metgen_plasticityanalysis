library(Hmsc)
library(colorspace)
library(vioplot)

localDir = "~/Documents/GitHub/Metgen_plasticityanalysis/Hmsc"
ModelDir = file.path(localDir, "models")
PanelDir = file.path(localDir, "panels")

TH=2
samples_list = c(10,250)
thin_list = c(1,10)
nst = length(thin_list)
nChains = 4
modelnames=c("CR", "AS")

for(i in 1:length(modelnames)){
  ma = NULL
  na = NULL
  ma_gamma = NULL
  na_gamma = NULL
  ma_rho = NULL
  na_rho = NULL
  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    filename = paste("model_",modelnames[i],"_TH",as.character(TH),"_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),'.Rdata',sep = "")
    load(file=file.path(ModelDir,filename))
    mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    psrf.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
    psrf.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      na = c(na,paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples)))
    }
    if(is.null(ma_gamma)){
      ma_gamma=psrf.gamma[,1]
      na_gamma = paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples))
    } else {
      ma_gamma = cbind(ma_gamma,psrf.gamma[,1])
      na_gamma = c(na_gamma,paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples)))
    }
    if(is.null(ma_gamma)){
      ma_rho=psrf.rho[,1]
      na_rho = paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples))
    } else {
      ma_rho = cbind(ma_rho,psrf.rho[,1])
      na_rho = c(na_rho,paste0("model_",as.character(i),"_",as.character(thin),",",as.character(samples)))
    }
  }
  panel.name=paste("MCMC_convergence_TH",as.character(TH),"_model_",as.character(i),"_lognormal.pdf",sep = "")
  pdf(file=file.path(PanelDir,panel.name))
  par(mfrow=c(2,1))
  vioplot(ma,col=rainbow_hcl(1),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
  vioplot(ma,col=rainbow_hcl(1),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
  par(mfrow=c(2,1))
  vioplot(ma_gamma,col=rainbow_hcl(1),names=na_gamma,ylim=c(0,max(ma_gamma)),main="psrf(gamma)")
  vioplot(ma_gamma,col=rainbow_hcl(1),names=na_gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
  par(mfrow=c(2,1))
  boxplot(ma_rho,col=rainbow_hcl(1),names=na_rho,ylim=c(0,max(ma_rho)),main="psrf(rho)")
  dev.off()
}

####Evaluate model fit

i=1
filename = paste("model_",modelnames[i],"_TH",as.character(TH),"_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),"_lognormal.Rdata",sep = "")
load(file=file.path(ModelDir,filename))
m_NoInt=m
m_NoInt_WAIC=computeWAIC(m_NoInt,byColumn =TRUE)
median(m_NoInt_WAIC,na.rm = TRUE)
predY.MF.NoInt=computePredictedValues(m_NoInt, expected=FALSE)
EpredY.MF.NoInt=apply(predY.MF.NoInt, MARGIN=1:2, mean)
par(mfrow=c(1))
plot(m_NoInt$Y,EpredY.MF.NoInt)
abline(coef = c(0,1),col="red")
MF.NoInt=evaluateModelFit(hM=m_NoInt, predY=predY.MF.NoInt)
mean(MF.NoInt$R2)

i=2
filename = paste("model_",modelnames[i],"_TH",as.character(TH),"_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),"_lognormal.Rdata",sep = "")
load(file=file.path(ModelDir,filename))
m_NoInt=m
m_NoInt_WAIC=computeWAIC(m_NoInt,byColumn =TRUE)
median(m_NoInt_WAIC,na.rm = TRUE)
predY.MF.NoInt=computePredictedValues(m_NoInt, expected=FALSE)
EpredY.MF.NoInt=apply(predY.MF.NoInt, MARGIN=1:2, mean)
par(mfrow=c(1))
plot(m_NoInt$Y,EpredY.MF.NoInt)
abline(coef = c(0,1),col="red")
MF.NoInt=evaluateModelFit(hM=m_NoInt, predY=predY.MF.NoInt)
mean(MF.NoInt$R2)
