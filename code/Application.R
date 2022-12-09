rm(list=ls())
#Step 1: load packages and functions--------------------------------------------
setwd("C:/Users/Lenovo/Desktop/git/STAT4011-Project-2")
PackageList =c('abind')
NewPackages=PackageList[!(PackageList %in%
                            installed.packages()[,"Package"])]
if(length(NewPackages)) install.packages(NewPackages)
lapply(PackageList,require,character.only=TRUE)

source("code/DGP.r")
source('code/EM_gaus.r')
source('code/EM_pois.R')
source("code/viterbi_gaus.R")
source("code/viterbi_pois.R")

#Step 2: Explore the Energy consumption data------------------------------------
samplesize=60
data4=read.csv("data/energy_rawdata.csv")
data4=data4[1:samplesize,]
x4=data4[,2]

#Step 2.1 model selection
gaus_mean_x4=mean(x4)
gaus_sd_last_x4=sqrt(var(x4))
len_table=length(2:4)+1
aic_bic_x4=array(NA,dim = c(2,len_table),dimnames = list(c('aic','bic'),c(paste0('m=',2:4),"selected m")))

for (m in 2:4){
  #m=3
  initP_last=rep(1/m,m)
  transP_last=matrix(0.3/m,m,m)+matrix(c(rep(0.7,m),rep(0,m*(m-1))),m,m)
  gaus_mean_last=rep(gaus_mean_x4,m)
  gaus_sd_last=rep(gaus_sd_last_x4,m)
  
  list_mstep_para_gaus=EM_gaus(T=length(x4),m,x4,gaus_mean_last,gaus_sd_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-10))
  
  aic_x4_m=list_mstep_para_gaus$aic
  bic_x4_m=list_mstep_para_gaus$bic
  aic_bic_x4[1,m-1]=aic_x4_m
  aic_bic_x4[2,m-1]=bic_x4_m
}

aic_bic_x4[1,len_table]=as.numeric(which.min(aic_bic_x4[1,])+1)
aic_bic_x4[2,len_table]=as.numeric(which.min(aic_bic_x4[2,])+1)
print(aic_bic_x4)

#Step 2.2 parameter estimation and State Classification
m=as.numeric(aic_bic_x4[1,len_table])
initP_last=rep(1/m,m)
transP_last=matrix(0.3/m,m,m)+matrix(c(rep(0.7,m),rep(0,m*(m-1))),m,m)
gaus_mean_last=rep(gaus_mean_x4,m)
gaus_sd_last=rep(gaus_sd_last_x4,m)

list_mstep_para_gaus=EM_gaus(T=length(x4),m,x4,gaus_mean_last,gaus_sd_last,
                             transP_last,initP_last,
                             num_ite=10^4,tol=10^(-3))

list_mstep_para_gaus
initP_est_x4=list_mstep_para_gaus$initP
transP_est_x4=list_mstep_para_gaus$transP
gaus_mean_est_x4=list_mstep_para_gaus$gaus_mean
gaus_var_est_x4=list_mstep_para_gaus$gaus_var

path=viterbi_gaus(initP_est_x4,transP_est_x4,gaus_mean_est_x4,gaus_var_est_x4,m,T=length(x4),x4)
path
#gaus_mean=gaus_mean_est_x2[order(gaus_mean_est_x2)]
path[path==1]=gaus_mean_est_x4[1]
path[path==2]=gaus_mean_est_x4[2]
path[path==3]=gaus_mean_est_x4[3]
#path[path==4]=gaus_mean[4]

#Plot the path
png("output/viterbi_energy.png",width=800,height=400)
par(mar = c(3, 3, 3, 8), xpd = TRUE)
plot(x4,type = "l",xlab = "time",ylab = "")
par(new=TRUE)
plot(path,type = "l",col="red",ylab = "state/data",yaxt='n',xlab = "time",
     main="estimated hidden state")
dev.off()


#Step 2.3 density compare-------------------------------------------------------
T=500
s=rep(NA,T)
xr=rep(NA,T)
p_t=initP_est_x4
s[1]=which.max(p_t)
xr[1]=rnorm(1,gaus_mean_est_x4[s[1]],sqrt(gaus_var_est_x4[s[1]]))
for (t in 2:T){
  p_t=initP_est_x4%*%transP_est_x4
  s[t]=which.max(p_t)
  xr[t]=rnorm(1,gaus_mean_est_x4[s[t]],sqrt(gaus_var_est_x4[s[t]]))
}

png("output/energy_density.png",width=800,height=400)
par(mar = c(2, 2, 2, 7), xpd = TRUE)
plot(density(x4))
par(new=TRUE)
plot(density(xr),col="red",yaxt="n",xaxt="n",
     main="comparison between simulation and real data")
legend("topright", inset=c(-0.25, 0),legend=c("real data","simulated data"),
       col=c("black", "red"),lty=1:2,cex = 0.6)
dev.off()


#Step 3: Explore the accident data----------------------------------------------
#Step 3.1 model select
samplesize=60
data3<-read.csv("data/accident_rawdata.csv")
x3<-data3[1:samplesize,]
c1<-as.numeric(gsub(",", "", x3$Total_Accident))
data3<-matrix(c1,ncol=1)
colnames(data3)="x"
x3=as.numeric(data3)

mean_x3=mean(x3)
len_table=length(2:4)+1
aic_bic_x3=array(NA,dim = c(2,len_table),dimnames = list(c('aic','bic'),c(paste0('m=',2:4),"selected m")))

for (m in 2:4){
  #m=2
  initP_last=rep(1/m,m)
  transP_last=matrix(0.3/m,m,m)+matrix(c(rep(0.7,m),rep(0,m*(m-1))),m,m)
  lambda_last=rep(mean_x3,m)
  
  list_mstep_para_pois=EM_pois(T=length(x3),m,x3,lambda_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-30))
  
  aic_x3_m=list_mstep_para_pois$aic
  bic_x3_m=list_mstep_para_pois$bic
  aic_bic_x3[1,m-1]=aic_x3_m
  aic_bic_x3[2,m-1]=bic_x3_m
}

aic_bic_x3[1,len_table]=as.numeric(which.min(aic_bic_x3[1,])+1)
aic_bic_x3[2,len_table]=as.numeric(which.min(aic_bic_x3[2,])+1)
print(aic_bic_x3)

#Step 3.2 parameter estimation and viterbi algorithm
m=4
initP_last=rep(1/m,m)
#transP_last=matrix(c(0.1,0.9,0.8,0.2),2,2,byrow = TRUE)
transP_last=matrix(0.3/m,m,m)+matrix(c(rep(0.7,m),rep(0,m*(m-1))),m,m)
lambda_last=rep(300,m)

list_mstep_para_pois=EM_pois(T=length(x3),m,x3,lambda_last,
                             transP_last,initP_last,
                             num_ite=10^4,tol=10^(-30))

lambda_est_x3=list_mstep_para_pois$lambda
transP_est_x3=list_mstep_para_pois$transP
initP_est_x3=list_mstep_para_pois$initP

path=viterbi_pois(initP_est_x3,transP_est_x3,lambda_est_x3,m,T=length(x3),x3)
path
#lambda_x1=lambda_est_x1[order(gaus_mean_est_x2)]
path[path==1]=lambda_est_x3[1]
path[path==2]=lambda_est_x3[2]
path[path==3]=lambda_est_x3[3]
#path[path==4]=gaus_mean[4]
path

png("output/viterbi_accident.png",width=800,height=400)
par(mar = c(3, 3, 3, 8), xpd = TRUE)
plot(x3,type = "l",xlab = "time",ylab = "")
par(new=TRUE)
plot(path,type = "l",col="red",ylab = "",yaxt='n',xlab = "")
legend("topright", inset=c(-0.4, 0),legend=c("data", "viterbi result"),
       col=c("black", "red"),lty=1:2,cex = 0.6)
dev.off()

#Step 3.3 density compare
T=500
s=rep(NA,T)
xr=rep(NA,T)
p_t=initP_est_x3
s[1]=which.max(p_t)
xr[1]=rpois(1,lambda_est_x3[s[1]])
for (t in 2:T){
  p_t=initP_est_x3%*%transP_est_x3
  s[t]=which.max(p_t)
  xr[t]=rpois(1,lambda_est_x3[s[t]])
}

png("output/accident_density.png",width=800,height=400)
par(mar = c(2, 2, 2, 7), xpd = TRUE)
plot(density(x3))
par(new=TRUE)
plot(density(xr),col="red",yaxt="n",xaxt="n")
legend("topright", inset=c(-0.25, 0),legend=c("original data", "estimated"),
       col=c("black", "red"),lty=1:2,cex = 0.6)
dev.off()