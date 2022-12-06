#Install packages and load functions-----------------------------
PackageList =c('abind')
NewPackages=PackageList[!(PackageList %in%
                            installed.packages()[,"Package"])]
if(length(NewPackages)) install.packages(NewPackages)
lapply(PackageList,require,character.only=TRUE)
setwd("C:/Users/Lenovo/Desktop/The Chinese University of HongKong/Year4T1/STAT4011/project 2")
# setwd("/Users/jiangyunhui/Downloads/STAT4011proj2")
source("DGP.r")
source('EM_gaus.r')
source('EM_pois.R')
source("viterbi_gaus.R")

#Simulation for Gaussian------------------------------------------
#Set up-----------------------------------------------------------
m=2
T=60
nRep=1000
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

initP_last=initP_true[1:(m-1)]+sample(noise_seq,size=(m-1))
initP_last=c(initP_last,1-sum(initP_last))
transP_last=transP_true+matrix(sample(noise_seq,size=m^2),
                                  nrow=m,ncol=m)
transP_last[1:m,m]=1-rowSums(transP_last[,-m])
gaus_sd_last=gaus_sd_true+sample(noise_seq,size=m)
gaus_mean_last=gaus_mean_true+sample(noise_seq,size=m)


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

#Conduct simulation through replications-------------------------
for(iRep in 1:nRep){
  set.seed(iRep)
  data=DGP_gaus(T,transP_true,initP_true,gaus_mean_true,gaus_sd_true)
  x=data$X
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
gaus_out=gaus_out[-which(is.na(gaus_out[,1])),]

#Examine MSE of estimators------------------------------------
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

plot(density(decoding[,T+1]),xlab="accuracy",ylab="density",
     main="density of accuracy over replications")