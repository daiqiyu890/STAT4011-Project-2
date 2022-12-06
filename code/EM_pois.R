#This script conduct E-M algorithm for Poisson distribution
# generate diags of P(x_t) for all x_t---------------------
dependP_diags_pois=function(T,m,x,lambda_last){
  xc_diag_m=c() # combine of row of P(x_t), jth row for diags of P(x_j)
  for (j in 1:T){
    xc_diag_r=c() # diag entries for P(x_t), ith entry is P(x_t|c_i)
    for (i in 1:m){
      xc_diag_r=c(xc_diag_r,dpois(x[j],lambda_last[i]))
    }
    xc_diag_m=rbind(xc_diag_m,xc_diag_r)
  }
  
  return(xc_diag_m)
}


# calculate forward-backward probability-----------------
Forward_Backward_pois=function(T,transP,initP,xc_diag_m,x,
                               lambda_last){
  
  # forward prob
  m=length(initP)
  Forward_m=function(xc_diag_m,T,initP,transP){
    alpha_t_ratio=rep(NA,T) # store the |alpha_t|/|alpha_t-1|, but 1st is |alpha_1|
    alpha_m=t(initP)%*%diag(xc_diag_m[1,]) # T*m matrix, i_th row stores alpha_t for t=i
    #print("alpha_m")
    #print(alpha_m)
    for (t in 2:T){
      alpha_t=t(initP)%*%diag(xc_diag_m[1,]) # vector of length m
      alpha_t_ratio[1]=sum(alpha_t)
      alpha_t=alpha_t/sum(alpha_t)
      for (tt in 2:t){
        alpha_t=alpha_t%*%transP%*%diag(xc_diag_m[tt,],m,m)
        alpha_t_ratio[tt]=sum(alpha_t)
        alpha_t=alpha_t/sum(alpha_t)
      }
      alpha_m=rbind(alpha_m,alpha_t)
    }
    list_alpha_m_ratio=list(alpha_m,alpha_t_ratio)
    names(list_alpha_m_ratio)=c("alpha_m","alpha_t_ratio")
    #print(list_alpha_m_ratio)
    return(list_alpha_m_ratio)
  }
  list_alpha_m_ratio=Forward_m(xc_diag_m,T,initP,transP)
  
  # backward prob
  Backward_m=function(xc_diag_m,transP){
    beta_m=matrix(NA,nrow=T,ncol=m) # T*m, i_th row stores beta_t for t=i
    for (t in 1:(T-1)){
      beta_t=diag(1,m)
      for (tt in (t+1):T){
        beta_t=beta_t%*%transP%*%diag(xc_diag_m[tt,],m)
        beta_t=beta_t/sum(beta_t)
      }
      beta_t=beta_t%*%(matrix(1,nrow=m,ncol=1))
      beta_m[t,]=c(beta_t)
    }
    #Tth beta
    beta_m[T,]=rep(1,m)
    list_beta_m=list(beta_m)
    names(list_beta_m)=c("beta_m")
    return(list_beta_m)
  }
  list_beta_m=Backward_m(xc_diag_m,transP)
  
  list_alpha_beta_m=append(list_alpha_m_ratio,list_beta_m)
  #print(list_alpha_beta_m)
  return(list_alpha_beta_m)
}

#E step----------------------------------------------------
E_step_pois=function(alpha_m,alpha_t_ratio,beta_m,T,transP,xc_diag_m){
  m=dim(transP)[1]
  L_T_m=matrix(NA,T,m)
  for (t in 1:T){
    L_T_m[t,]=rep(as.numeric(t(alpha_m[t,])%*%beta_m[t,]),m)
  }
  u_jt_m=alpha_m*beta_m/L_T_m #T*m, i_th row stores hat(u_j(t)) for t=i
  # v_jkt
  v_jkt_m=c() #3d array for hat(v_jk(t)), jkt:123 (alpha,beta,t), t: 2:T
  for (t in 2:T){
    v_jk_t=(alpha_m[t-1,]%*%t(beta_m[t,]))*transP #jk:alpha_t-1(j)*beta_t+1(k)
    p_kt_m=matrix(rep(xc_diag_m[t,],m),m,m,byrow=TRUE)
    v_jk_t=v_jk_t*p_kt_m/(as.numeric(t(alpha_m[t,])%*%beta_m[t,]))/alpha_t_ratio[t]
    v_jkt_m=abind(v_jkt_m,v_jk_t,along=3)
  }
  list_uv=list(u_jt_m,v_jkt_m)
  names(list_uv)=c("u_jt_m","v_jkt_m")
  return(list_uv)
}
#M step-----------------------------------------------------
M_step_pois=function(u_jt_m,v_jkt_m,x,m,T){
  # initial distribution
  initP=u_jt_m[1,]/sum(u_jt_m[1,])
  
  # transition prob
  f_jk=apply(v_jkt_m,c(1,2),sum)
  #print("f_jk")
  #print(f_jk)
  denominator_m=t(matrix(rep(rowSums(f_jk,dims = 1),m),m,m,byrow = TRUE))
  transP=f_jk/denominator_m
  
  # dependP
  x_m=t(matrix(rep(x,m),m,T,byrow=TRUE))
  
  # pois_lambda
  lambda=colSums(x_m*u_jt_m)/colSums(u_jt_m)
  list_mstep_para_pois=list(initP,transP,lambda)
  
  #update them according to certain order
  order_pois=order(lambda)
  lambda=lambda[order_pois]
  initP=initP[order_pois]
  transP_new=transP
  for(i in 1:m){
    for(j in 1:m){
      transP_new[order_pois[i],order_pois[j]]=transP[i,j]
    }
  }
  transP=transP_new
  
  names(list_mstep_para_pois)=c("initP","transP","lambda")
  return(list_mstep_para_pois)
}


EM_pois=function(T,m,x,lambda_last,
                 transP_last,initP_last,num_ite=10^3,tol=10^(-3)){
  loglikelihood=rep(NA,num_ite)
  for(i_ite in 1:num_ite){
    xc_diag_m=
      dependP_diags_pois(T,m,x,lambda_last)
    
    list_alpha_beta_m=
      Forward_Backward_pois(T,transP_last,initP_last,xc_diag_m,x,
                            lambda_last)
    
    alpha_m=list_alpha_beta_m$alpha_m
    alpha_t_ratio=list_alpha_beta_m$alpha_t_ratio
    beta_m=list_alpha_beta_m$beta_m
    
    list_uv=E_step_pois(alpha_m,alpha_t_ratio,beta_m,T,transP_last,xc_diag_m)
    u_jt_m=list_uv$u_jt_m
    v_jkt_m=list_uv$v_jkt_m
    
    list_mstep_para_pois=M_step_pois(u_jt_m,v_jkt_m,x,m,T)
    
    #Update estimates
    lambda_last=list_mstep_para_pois$lambda
    transP_last=list_mstep_para_pois$transP
    initP_last=list_mstep_para_pois$initP
  
    #Likelihood
    loglikelihood[i_ite]=log(as.numeric(t(alpha_m[T,])%*%beta_m[T,]))+sum(log(alpha_t_ratio))
    if(is.na(loglikelihood[i_ite])){
      break
    }
    if(i_ite>1){
      if(abs(loglikelihood[i_ite]-loglikelihood[i_ite-1])<tol){
        break
      }
    }
  }
  if(is.na(loglikelihood[i_ite])){
    list_mstep_para_gaus=NULL
    return(list_mstep_para_pois)
  } else {
    list_mstep_para_pois$loglikelihood=loglikelihood[i_ite]
    p=m^2+m-1
    list_mstep_para_pois$aic=-2*list_mstep_para_pois$loglikelihood+2*p
    list_mstep_para_pois$bic=-2*list_mstep_para_pois$loglikelihood+p*log(T)
    return(list_mstep_para_pois) 
  }
}