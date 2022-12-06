setwd("C:/Users/Lenovo/Desktop/git/STAT4011-Project-2")
#Step 1: data pre-processing----------------------------------------------------
#visitors
data1<-read.csv("data/visitors_rawdata.csv")
samplesize=60
a<-nrow(data1)
b<-a-(samplesize-1)
x1<-data1[b:a,]

c1<-as.numeric(gsub(",", "", x1$Unique.Visits))
Unique.Visits<-c1

x1<-as.numeric(gsub(",", "", x1$Unique.Visits))
mean1<-mean(c1)
sd1<-sd(c1)
Unique.Visits<-(c1-mean1)/sd1

data1<-matrix(Unique.Visits,ncol=1)
colnames(data1)="x"

#sales
data2<-read.csv("data/sales_rawdata.csv")
a<-nrow(data2)
b<-a-(samplesize-1)
x2<-data2[b:a,]
x2=x2$ProductP1
x2=(x2-mean(x2))/sd(x2)


#Step 2: Model Selection--------------------------------------------------------
source("code/DGP.r")
source('code/EM_gaus.r')
source('code/EM_pois.R')
source("code/viterbi_gaus.R")
source("code/viterbi_pois.R")

#sales

x2=c(data2$ProductP4)
gaus_mean_x2=mean(x2)
gaus_sd_last_x2=sqrt(var(x2))
len_table=length(2:4)+1
aic_bic_x2=array(NA,dim = c(2,len_table),dimnames = list(c('aic','bic'),c(paste0('m=',2:4),"selected m")))

for (m in 2:4){
  #m=3
  initP_last=rep(1/m,m)
  transP_last=matrix(1/m,m,m)
  gaus_mean_last=rep(gaus_mean_x2,m)
  gaus_sd_last=rep(gaus_sd_last_x2,m)
  
  list_mstep_para_gaus=EM_gaus(T=length(x2),m,x2,gaus_mean_last,gaus_sd_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-3))

  aic_x2_m=list_mstep_para_gaus$aic
  bic_x2_m=list_mstep_para_gaus$bic
  aic_bic_x2[1,m-1]=aic_x2_m
  aic_bic_x2[2,m-1]=bic_x2_m
}

aic_bic_x2[1,len_table]=as.numeric(which.min(aic_bic_x2[1,])+1)
aic_bic_x2[2,len_table]=as.numeric(which.min(aic_bic_x2[2,])+1)
print(aic_bic_x2)
#pois==============================
data1<-read.csv("/Users/jiangyunhui/Documents/GitHub/STAT4011-Project-2/data/visitors_rawdata.csv")
samplesize=60
a<-nrow(data1)
b<-a-(samplesize-1)
x1<-data1[b:a,]
x1<-as.numeric(gsub(",", "", x1$Unique.Visits))
mean_x1=mean(x1)
len_table=length(2:4)+1
aic_bic_x1=array(NA,dim = c(2,len_table),dimnames = list(c('aic','bic'),c(paste0('m=',2:4),"selected m")))

for (m in 2:4){
  #m=3
  initP_last=rep(1/m,m)
  transP_last=matrix(1/m,m,m)
  lambda_x1=rep(mean_x1,m)
  
  list_mstep_para_pois=EM_pois(T=length(x1),m,x1,lambda_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-30))
  
  aic_x1_m=list_mstep_para_pois$aic
  bic_x1_m=list_mstep_para_pois$bic
  aic_bic_x1[1,m-1]=aic_x1_m
  aic_bic_x1[2,m-1]=bic_x1_m
}

aic_bic_x1[1,len_table]=as.numeric(which.min(aic_bic_x1[1,])+1)
aic_bic_x1[2,len_table]=as.numeric(which.min(aic_bic_x1[2,])+1)
print(aic_bic_x1)

m=2
initP_last=rep(1/m,m)
transP_last=matrix(c(0.1,0.9,0.8,0.2),2,2)
lambda_last=rep(3000,m)

list_mstep_para_pois=EM_pois(T=length(x1),m,x1,lambda_last,
                             transP_last,initP_last,
                             num_ite=10^4,tol=10^(-30))

lambda_est_x1=list_mstep_para_pois$lambda
transP_est_x1=list_mstep_para_pois$transP
initP_est_x1=list_mstep_para_pois$initP

path=viterbi_pois(initP_est_x1,transP_est_x1,lambda_est_x1,m,T=length(x1),x1)
path
#lambda_x1=lambda_est_x1[order(gaus_mean_est_x2)]
path[path==1]=lambda_est_x1[1]
path[path==2]=lambda_est_x1[2]
#path[path==3]=gaus_mean[3]
#path[path==4]=gaus_mean[4]
path
plot(x1,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

#accident===============================

data3<-read.csv("/Users/jiangyunhui/Documents/GitHub/STAT4011-Project-2/data/accident.csv")
samplesize=60
x3=data3$x
mean_x3=mean(x3)
len_table=length(2:4)+1
aic_bic_x3=array(NA,dim = c(2,len_table),dimnames = list(c('aic','bic'),c(paste0('m=',2:4),"selected m")))

for (m in 2:4){
  #m=2
  initP_last=rep(1/m,m)
  transP_last=matrix(1/m,m,m)
  lambda_x3=rep(mean_x3,m)
  
  list_mstep_para_pois=EM_pois(T=length(x3),m,x3,lambda_last,
                               transP_last,initP_last,
                               num_ite=10^4,tol=10^(-30))
  
  aic_x3_m=list_mstep_para_pois$aic
  bic_x3_m=list_mstep_para_pois$bic
  aic_bic_x3[1,m-1]=aic_x1_m
  aic_bic_x3[2,m-1]=bic_x1_m
}

aic_bic_x3[1,len_table]=as.numeric(which.min(aic_bic_x3[1,])+1)
aic_bic_x3[2,len_table]=as.numeric(which.min(aic_bic_x3[2,])+1)
print(aic_bic_x3)

m=3
initP_last=rep(1/m,m)
#transP_last=matrix(c(0.1,0.9,0.8,0.2),2,2,byrow = TRUE)
transP_last=matrix(c(0.1,0.7,0.2,0.8,0.1,0.1,0.1,0.6,0.3),3,3,byrow = TRUE)
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
plot(x3,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

#Step 3: Parameter Estimation---------------------------------------------------

m=as.numeric(aic_bic_x2[1,len_table])

#m=as.numeric(aic_bic_x2[1,len_table])
m=4

initP_last=rep(1/m,m)
transP_last=matrix(1/m,m,m)
gaus_mean_last=rep(10,m)
gaus_sd_last=rep(10,m)

list_mstep_para_gaus_est=EM_gaus(T=length(x2),m,x2,gaus_mean_last,gaus_sd_last,
                             transP_last,initP_last,
                             num_ite=10^4,tol=10^(-200))
                             num_ite=10^4,tol=10^(-3))

gaus_mean_est_x2=list_mstep_para_gaus_est$gaus_mean
gaus_var_est_x2=list_mstep_para_gaus_est$gaus_var
transP_est_x2=list_mstep_para_gaus_est$transP
initP_est_x2=list_mstep_para_gaus_est$initP

# viterbi & plot
path=viterbi_gaus(initP_est_x2,transP_est_x2,gaus_mean_est_x2,gaus_var_est_x2,m,T=length(x2),x2)
path
gaus_mean=gaus_mean_est_x2[order(gaus_mean_est_x2)]
path[path==1]=gaus_mean[1]
path[path==2]=gaus_mean[2]
#path[path==3]=gaus_mean[3]
#path[path==4]=gaus_mean[4]
path
plot(x2,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

#Step 4: Forecast---------------------------------------------------------------


#Step 4: Plot-------------------------------------------------------------------





