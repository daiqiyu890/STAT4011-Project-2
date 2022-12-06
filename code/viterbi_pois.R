# viterbi

viterbi_pois=function(initP,transP,lambda,m,T,x){
  v=matrix(rep(NA,T*m),T,m)
  back_pointer=matrix(rep(NA,T*m),T,m)
  
  dependP_x_m=matrix(rep(NA,T*m),T,m)
  for (t in 1:T){
    for (j in 1:m){
      dependP_x_m[t,j]=dpois(x[t],lambda[j])
    }
  }
  #initialize
  v[1,]=t(initP)*t(dependP_x_m[1,])
  back_pointer[1,]=rep(0,m)
  #recursion
  for (t in 2:T){
    for (j in 1:m){
      ker=v[t-1,]*transP[,j]
      back_pointer[t,j]=which.max(ker)
      v[t,j]=ker[back_pointer[t,j]]*dependP_x_m[t,j]
    }
    v[t,]=v[t,]/sum(v[t,])
  }
  #termination
  best_path=rep(NA,T)
  best_path[T]=which.max(v[T,])
  for (t in 1:(T-1)){
    best_path[T-t]=back_pointer[T-t+1,best_path[T-t+1]]
  }
  
  #list_bestpath_bestpathpro=list(best_path,)
  #names(list_bestpath_bestpathpro)=c("best_path","best_path_pro")
  return(best_path)
}

