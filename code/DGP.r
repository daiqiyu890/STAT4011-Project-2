DGP_gaus=function(T,transP,initP,gaus_mean,gaus_sd){
  m=length(initP)
  C=X=numeric(T)
  C[1]=sample(1:m,size=1,prob=initP)
  X[1]=rnorm(n=1,mean=gaus_mean[C[1]],sd=gaus_sd[C[1]])
  for(it in 2:T){
    p=transP[C[it-1],]
    C[it]=sample(1:m,size=1,prob=p)
    X[it]=rnorm(n=1,mean=gaus_mean[C[it]],sd=gaus_sd[C[it]])
  }
  out=list(C,X)
  names(out)=c("C","X")
  return(out)
}

DGP_pois=function(T,transP,initP,lambda_pois){
  m=length(initP)
  C=X=numeric(T)
  C[1]=sample(1:m,size=1,prob=initP)
  X[1]=rpois(n=1,lambda=lambda_pois[C[1]])
  for(it in 2:T){
    p=transP[C[it-1],]
    C[it]=sample(1:m,size=1,prob=p)
    X[it]=rpois(n=1,lambda=lambda_pois[C[it]])
  }
  out=list(C,X)
  names(out)=c("C","X")
  return(out)
}

