rm(list=ls())
#Step 1: Install packages and load functions------------------------------------
PackageList =c('abind')
NewPackages=PackageList[!(PackageList %in%
                            installed.packages()[,"Package"])]
if(length(NewPackages)) install.packages(NewPackages)
lapply(PackageList,require,character.only=TRUE)
setwd("C:/Users/Lenovo/Desktop/git/STAT4011-Project-2")
# setwd("/Users/jiangyunhui/Downloads/STAT4011proj2")
source("code/DGP.r")
source('code/EM_gaus.r')
source('code/EM_pois.R')
source("code/viterbi_gaus.R")
source("code/viterbi_pois.R")


#Step 2: Simulation for Gaussian------------------------------------------------
#Step 2.1: Set up
m=4
T=60
nRep=10^3
#ture parameters
#Note that we require the mean of different states with increasing order
if(m==2){
  initP_true=c(0.2,0.8) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.8,0.2,0.25,0.75),nrow=m,byrow=TRUE)
  gaus_sd_true=sqrt(c(1,4))
  gaus_mean_true=c(0,4)
} else if(m==3){
  initP_true=c(0.1,0.3,0.6) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.1,0.2,0.7,
                       0.25,0.5,0.25,
                       0.4,0.15,0.45),nrow=m,byrow=TRUE)
  gaus_sd_true=sqrt(c(1,4,2))
  gaus_mean_true=c(0,2,4)  
} else if(m==4){
  initP_true=c(0.1,0.3,0.2,0.4) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.1,0.2,0.25,0.45,
                       0.15,0.4,0.15,0.3,
                       0.22,0.32,0.2,0.26,
                       0.08,0.42,0.3,0.2),
                     nrow=m,byrow=TRUE)
  gaus_sd_true=sqrt(c(1,4,2,3))
  gaus_mean_true=c(0,2,3,4)  
}

gaus_para_true=c(gaus_mean_true,gaus_sd_true^2,initP_true,
                 as.vector(transP_true))

#initial parameters
noise_seq=seq(0.001,0.01,length=100)
set.seed(4011)

transP_name=c()
for(i in 1:m){
  transP_name=c(transP_name,paste0("transP",1:m,i))
}
gaus_out=array(NA,dim=c(nRep,3+m*3+m^2),
               dimnames=list(paste0("iRep=",1:nRep),
                             c("aic","bic","llk",
                               paste0("mean",1:m),
                               paste0("var",1:m),paste0("init",1:m),
                               transP_name)))

decoding=array(NA,dim=c(nRep,T+1),
               dimnames=list(paste0("iRep=",1:nRep),
                             c(paste0("t=",1:T),"accuracy")))

#Step 2.2 Conduct simulation through replications
for(iRep in 1:nRep){
  set.seed(iRep)
  data=DGP_gaus(T,transP_true,initP_true,gaus_mean_true,gaus_sd_true)
  x=data$X
  
  initP_last=initP_true[1:(m-1)]+sample(noise_seq,size=(m-1))
  initP_last=c(initP_last,1-sum(initP_last))
  transP_last=transP_true+matrix(sample(noise_seq,size=m^2),
                                 nrow=m,ncol=m)
  transP_last[1:m,m]=1-rowSums(as.matrix(transP_last[,-m],nrow=m))
  
  gaus_sd_last=gaus_sd_true+sample(noise_seq,size=m)
  gaus_mean_last=gaus_mean_true+sample(noise_seq,size=m)
  list_mstep_para_gaus=EM_gaus(T,m,x,gaus_mean_last,gaus_sd_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-3))
  if(!is.null(list_mstep_para_gaus)){
    gaus_out[iRep,"aic"]=list_mstep_para_gaus$aic
    gaus_out[iRep,"bic"]=list_mstep_para_gaus$bic
    gaus_out[iRep,"llk"]=list_mstep_para_gaus$loglikelihood
    gaus_out[iRep,4:(3+m)]=list_mstep_para_gaus$gaus_mean
    gaus_out[iRep,(4+m):(3+2*m)]=list_mstep_para_gaus$gaus_var
    gaus_out[iRep,(4+2*m):(3+3*m)]=list_mstep_para_gaus$initP
    gaus_out[iRep,(4+3*m):(3+3*m+m^2)]=as.vector(list_mstep_para_gaus$transP)
    decoding[iRep,1:T]=viterbi_gaus(list_mstep_para_gaus$initP,
                                    list_mstep_para_gaus$transP,
                                    list_mstep_para_gaus$gaus_mean,
                                    list_mstep_para_gaus$gaus_var,
                                    m,T,x)
    decoding[iRep,(T+1)]=sum(decoding[iRep,1:T]==data$C)/T
  }
  print(iRep)
}

if(length(which(is.na(gaus_out[,1])))>0){
  gaus_out=gaus_out[-which(is.na(gaus_out[,1])),]
}

if(length(which(is.na(decoding[,(T+1)])))>0){
  decoding=decoding[-which(is.na(decoding[,(T+1)])),]
}

#Step 2.3 Examine MSE of estimators
gaus_MSE=array(dim=c(4,m*3+m^2),dimnames=list(c("truth","bias","variance","MSE"),
                                      c(paste0("mean",1:m),
                                          paste0("var",1:m),paste0("init",1:m),
                                          transP_name)))
gaus_MSE[1,]=gaus_para_true
gaus_MSE[2,]=apply(gaus_out[,-(1:3)],MARGIN=2,FUN=mean)-gaus_para_true
gaus_MSE[3,]=apply(gaus_out[,-(1:3)],MARGIN=2,FUN=var)
gaus_MSE[4,]=(gaus_MSE[2,])^2+gaus_MSE[3,]
file_name=paste0("output/gaus_MSE_m",m,".csv")
write.csv(gaus_MSE,file_name)

file_name=paste0("output/gaus_accuracy_m",m,".csv")
write.csv(decoding,file_name)

file_name=paste0("output/gaus_accuracy_m",m,".png")
png(file_name)
plot(density(decoding[,T+1]),xlab="accuracy",ylab="density",
     main="density of accuracy over replications")
dev.off()

#Step 2.4 Test AIC BIC, in terms of model selection
model_sele=function(m_test){
  if(m_test==2){
    initP_true=c(0.2,0.8) #initial distribution of Hidden Markov States
    transP_true=matrix(c(0.8,0.2,0.25,0.75),nrow=m_test,byrow=TRUE)
    gaus_sd_true=sqrt(c(1,4))
    gaus_mean_true=c(0,4)
  } else if(m_test==3){
    initP_true=c(0.1,0.3,0.6) #initial distribution of Hidden Markov States
    transP_true=matrix(c(0.1,0.2,0.7,
                         0.25,0.5,0.25,
                         0.4,0.15,0.45),nrow=m_test,byrow=TRUE)
    gaus_sd_true=sqrt(c(1,4,2))
    gaus_mean_true=c(0,2,4)  
  } else if(m_test==4){
    initP_true=c(0.1,0.3,0.2,0.4) #initial distribution of Hidden Markov States
    transP_true=matrix(c(0.1,0.2,0.25,0.45,
                         0.15,0.4,0.15,0.3,
                         0.22,0.32,0.2,0.26,
                         0.08,0.42,0.3,0.2),
                       nrow=m_test,byrow=TRUE)
    gaus_sd_true=sqrt(c(1,4,2,3))
    gaus_mean_true=c(0,2,3,4)  
  }
  
  initP_last=initP_true[1:(m_test-1)]+sample(noise_seq,size=(m_test-1))
  initP_last=c(initP_last,1-sum(initP_last))
  transP_last=transP_true+matrix(sample(noise_seq,size=m_test^2),
                                 nrow=m_test,ncol=m_test)
  
  transP_last[1:m_test,m_test]=1-rowSums(as.matrix(transP_last[,-m_test],nrow=m_test))
  gaus_sd_last=gaus_sd_true+sample(noise_seq,size=m_test)
  gaus_mean_last=gaus_mean_true+sample(noise_seq,size=m_test)
  
  list_mstep_para_gaus=EM_gaus(T,m_test,x,gaus_mean_last,gaus_sd_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-3))
  return(list_mstep_para_gaus)
}

best_model=array(NA,dim=c(nRep,8),dimnames=list(paste0("iRep=",1:nRep),
                                               c("AIC_2","BIC_2","AIC_3","BIC_3",
                                                 "AIC_4","BIC_4","AIC_m","BIC_m")))


for(iRep in 1:nRep){
  set.seed(iRep)
  data=DGP_gaus(T,transP_true,initP_true,gaus_mean_true,gaus_sd_true)
  x=data$X
  
  list_mstep_para_gaus2=model_sele(m_test=2)
  list_mstep_para_gaus3=model_sele(m_test=3)
  list_mstep_para_gaus4=model_sele(m_test=4)
  if(!(is.null(list_mstep_para_gaus2) | is.null(list_mstep_para_gaus3) | 
      is.null(list_mstep_para_gaus4))){
    best_model[iRep,"AIC_2"]=list_mstep_para_gaus2$aic
    best_model[iRep,"BIC_2"]=list_mstep_para_gaus2$bic
    

    best_model[iRep,"AIC_3"]=list_mstep_para_gaus3$aic
    best_model[iRep,"BIC_3"]=list_mstep_para_gaus3$bic

    best_model[iRep,"AIC_4"]=list_mstep_para_gaus4$aic
    best_model[iRep,"BIC_4"]=list_mstep_para_gaus4$bic
    
    best_model[iRep,"AIC_m"]=as.numeric(which.min(best_model[iRep,c(1,3,5)])+1)
    best_model[iRep,"BIC_m"]=as.numeric(which.min(best_model[iRep,c(2,4,6)])+1)
  }
  
  print(iRep)
}

if(length(which(is.na(best_model[,1])))){
  best_model=best_model[-which(is.na(best_model[,1])),]
}

file_name=paste0("output/model_sele_m=",m,".csv")
write.csv(best_model,file=file_name)
print(paste0("AIC select the correct model in ",
             sum(best_model[,7]==m)*100/nrow(best_model),"% cases"))

print(paste0("BIC select the correct model in ",
             sum(best_model[,8]==m)*100/nrow(best_model),"% cases"))


#Step 3 Simulation for Poisson--------------------------------------------------
#Step 3.1 Set up
m=2
T=60
nRep=10^3
#ture parameters
#Note that we require the mean of different states with increasing order
if(m==2){
  initP_true=c(0.2,0.8) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.8,0.2,0.25,0.75),nrow=m,byrow=TRUE)
  lambda_true=c(1,3)
} else if(m==3){
  initP_true=c(0.1,0.3,0.6) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.1,0.2,0.7,
                       0.25,0.5,0.25,
                       0.4,0.15,0.45),nrow=m,byrow=TRUE)
  lambda_true=c(1,3,5)
} else if(m==4){
  initP_true=c(0.1,0.3,0.2,0.4) #initial distribution of Hidden Markov States
  transP_true=matrix(c(0.1,0.2,0.25,0.45,
                       0.15,0.4,0.15,0.3,
                       0.22,0.32,0.2,0.26,
                       0.08,0.42,0.3,0.2),
                     nrow=m,byrow=TRUE)
  lambda_true=c(1,3,5,7)
}

pois_para_true=c(lambda_true,initP_true,
                 as.vector(transP_true))

#initial parameters
noise_seq=seq(0.001,0.01,length=100)
set.seed(4011)

transP_name=c()
for(i in 1:m){
  transP_name=c(transP_name,paste0("transP",1:m,i))
}
pois_out=array(NA,dim=c(nRep,3+m*2+m^2),
               dimnames=list(paste0("iRep=",1:nRep),
                             c("aic","bic","llk",
                               paste0("lambda",1:m),paste0("init",1:m),
                               transP_name)))

decoding=array(NA,dim=c(nRep,T+1),
               dimnames=list(paste0("iRep=",1:nRep),
                             c(paste0("t=",1:T),"accuracy")))

#Step 3.2 Conduct simulation through replications
for(iRep in 1:nRep){
  set.seed(iRep)
  data=DGP_pois(T,transP_true,initP_true,lambda_true)
  x=data$X
  
  m=2
  initP_last=initP_true[1:(m-1)]+sample(noise_seq,size=(m-1))
  initP_last=c(initP_last,1-sum(initP_last))
  transP_last=transP_true+matrix(sample(noise_seq,size=m^2),
                                 nrow=m_test,ncol=m_test)

  transP_last[1:m,m]=1-rowSums(as.matrix(transP_last[,-m],nrow=m))
  lambda_last=lambda_true+sample(noise_seq,size=m)

  list_mstep_para_pois=EM_pois(T,m,x,lambda_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-3))
  if(!is.null(list_mstep_para_pois)){
    pois_out[iRep,"aic"]=list_mstep_para_pois$aic
    pois_out[iRep,"bic"]=list_mstep_para_pois$bic
    pois_out[iRep,"llk"]=list_mstep_para_pois$loglikelihood
    pois_out[iRep,4:(3+m)]=list_mstep_para_pois$lambda
    pois_out[iRep,(4+m):(3+2*m)]=list_mstep_para_pois$initP
    pois_out[iRep,(4+2*m):(3+2*m+m^2)]=as.vector(list_mstep_para_pois$transP)
 
    
    decoding[iRep,1:T]=viterbi_pois(list_mstep_para_pois$initP,
                                    list_mstep_para_pois$transP,
                                    list_mstep_para_pois$lambda,
                                    m,T,x)
    decoding[iRep,(T+1)]=sum(decoding[iRep,1:T]==data$C)/T
  }
  print(iRep)
}

if(length(which(is.na(pois_out[,1])))>0){
  pois_out=pois_out[-which(is.na(pois_out[,1])),]
}

if(length(which(is.na(decoding[,(T+1)])))>0){
  decoding=decoding[-which(is.na(decoding[,(T+1)])),]
}

#Step 3.3 Examine MSE of estimators
pois_MSE=array(dim=c(4,m*2+m^2),dimnames=list(c("truth","bias","variance","MSE"),
                                              c(paste0("lambda",1:m),
                                                paste0("init",1:m),
                                                transP_name)))
pois_MSE[1,]=pois_para_true
pois_MSE[2,]=apply(pois_out[,-(1:3)],MARGIN=2,FUN=mean)-pois_para_true
pois_MSE[3,]=apply(pois_out[,-(1:3)],MARGIN=2,FUN=var)
pois_MSE[4,]=(pois_MSE[2,])^2+pois_MSE[3,]
file_name=paste0("output/pois_MSE_m",m,".csv")
write.csv(pois_MSE,file_name)

file_name=paste0("output/pois_accuracy_m",m,".csv")
write.csv(decoding,file_name)

file_name=paste0("output/pois_accuracy_m",m,".png")
png(file_name)
plot(density(decoding[,T+1]),xlab="accuracy",ylab="density",
     main="density of accuracy over replications")
dev.off()
