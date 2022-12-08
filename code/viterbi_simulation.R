# viterbi
# pois
setwd("/Users/jiangyunhui/Documents/GitHub/STAT4011-Project-2")
source("DGP.r")
source('EM_gaus.r')
source('EM_pois.R')
source("viterbi_gaus.R")
source("viterbi_pois.R")

#===============================================================
#pois
T=400
transP_true=matrix(c(0.9,0.02,0.05,0.03,
                     0.1,0.75,0.06,0.09,
                     0.05,0.1,0.55,0.3,
                     0.1,0.1,0.65,0.15),4,4,byrow = TRUE)

initP_true=c(0.2,0.4,0.1,0.3)
m=length(initP_true)
lambda_true=c(5,10,15,20)
data=DGP_pois(T,transP_true,initP_true,lambda_true)
x=data$X

# m=4
initP_last=c(0.5,0.2,0.2,0.1)
transP_last=matrix(rep(0.25,16),4,4,byrow = TRUE)
lambda_last=c(20,15,20,15)

list_mstep_para_pois=EM_pois(T,m,x,lambda_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_pois

lambda=list_mstep_para_pois$lambda
transP=list_mstep_para_pois$transP
initP=list_mstep_para_pois$initP

path=viterbi_pois(initP,transP,lambda,m,T,x)
path
path[path==1]=lambda[1]
path[path==2]=lambda[2]
path[path==3]=lambda[3]
path[path==4]=lambda[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

c=data$C
plot(c,type="l")
par(new=TRUE)
plot(path,type = "l",col="red")

#m=3
m=3
initP_last=c(0.5,0.2,0.3)
transP_last=matrix(rep(1/3,9),3,3,byrow = TRUE)
lambda_last=c(20,15,20)

list_mstep_para_pois=EM_pois(T,m,x,lambda_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_pois

lambda=list_mstep_para_pois$lambda
transP=list_mstep_para_pois$transP
initP=list_mstep_para_pois$initP

path=viterbi_pois(initP,transP,lambda,m,T,x)
path
path[path==1]=lambda[1]
path[path==2]=lambda[2]
path[path==3]=lambda[3]
#path[path==4]=lambda[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

#m=2
m=2
initP_last=c(0.5,0.5)
transP_last=matrix(rep(1/2,4),2,2,byrow = TRUE)
lambda_last=c(20,15)

list_mstep_para_pois=EM_pois(T,m,x,lambda_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_pois

lambda=list_mstep_para_pois$lambda
transP=list_mstep_para_pois$transP
initP=list_mstep_para_pois$initP

path=viterbi_pois(initP,transP,lambda,m,T,x)
path
path[path==1]=lambda[1]
path[path==2]=lambda[2]
#path[path==3]=lambda[3]
#path[path==4]=lambda[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")
#==============================================================
#gaus
T=400
transP_true=matrix(c(0.9,0.02,0.05,0.03,
                     0.1,0.75,0.06,0.09,
                     0.05,0.1,0.55,0.3,
                     0.1,0.1,0.65,0.15),4,4,byrow = TRUE)

initP_true=c(0.2,0.4,0.1,0.3)
m=length(initP_true)
gaus_mean_true=c(3,8,15,20)
gaus_sd_true=c(5,10,15,20)
data=DGP_gaus(T,transP_true,initP_true,gaus_mean_true,gaus_sd_true)
x=data$X

# m=4
initP_last=c(0.5,0.2,0.2,0.1)
transP_last=matrix(rep(0.25,16),4,4,byrow = TRUE)
gaus_mean_last=c(20,15,20,15)
gaus_sd_last=rep(5,4)

list_mstep_para_gaus=EM_gaus(T,m,x,gaus_mean_last,gaus_sd_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_gaus

gaus_mean=list_mstep_para_gaus$gaus_mean
gaus_var=list_mstep_para_gaus$gaus_var
transP=list_mstep_para_gaus$transP
initP=list_mstep_para_gaus$initP

path=viterbi_gaus(initP,transP,gaus_mean,gaus_var,m,T,x)
path
gaus_mean=gaus_mean[order(gaus_mean)]
path[path==1]=gaus_mean[1]
path[path==2]=gaus_mean[2]
path[path==3]=gaus_mean[3]
path[path==4]=gaus_mean[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

c=data$C
plot(c,type="l")
par(new=TRUE)
plot(path,type = "l",col="red")

#m=3
m=3
initP_last=c(0.5,0.2,0.3)
transP_last=matrix(rep(1/3,9),3,3,byrow = TRUE)
gaus_mean_last=rep(10,3)
gaus_sd_last=rep(5,3)

list_mstep_para_gaus=EM_gaus(T,m,x,gaus_mean_last,gaus_sd_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_gaus

gaus_mean=list_mstep_para_gaus$gaus_mean
gaus_var=list_mstep_para_gaus$gaus_var
transP=list_mstep_para_gaus$transP
initP=list_mstep_para_gaus$initP

path=viterbi_gaus(initP,transP,gaus_mean,gaus_var,m,T,x)
path
gaus_mean=gaus_mean[order(gaus_mean)]
path[path==1]=gaus_mean[1]
path[path==2]=gaus_mean[2]
path[path==3]=gaus_mean[3]
#path[path==4]=gaus_mean[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")

#m=2
m=2
initP_last=c(0.5,0.5)
transP_last=matrix(rep(1/2,4),2,2,byrow = TRUE)
gaus_mean_last=rep(10,2)
gaus_sd_last=rep(5,2)

list_mstep_para_gaus=EM_gaus(T,m,x,gaus_mean_last,gaus_sd_last,
                             transP_last,initP_last,
                             num_ite=10^3,tol=10^(-5))
list_mstep_para_gaus

gaus_mean=list_mstep_para_gaus$gaus_mean
gaus_var=list_mstep_para_gaus$gaus_var
transP=list_mstep_para_gaus$transP
initP=list_mstep_para_gaus$initP

path=viterbi_gaus(initP,transP,gaus_mean,gaus_var,m,T,x)
path
gaus_mean=gaus_mean[order(gaus_mean)]
path[path==1]=gaus_mean[1]
path[path==2]=gaus_mean[2]
#path[path==3]=gaus_mean[3]
#path[path==4]=gaus_mean[4]
path
plot(x,type = "l")
par(new=TRUE)
plot(path,type = "l",col="red")
