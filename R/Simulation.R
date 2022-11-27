X_gen=function(cov,x){
  chol_cov=chol(cov)
  return (x%*%chol_cov)
}

Cov=function(r,c,min,max){
  x=runif(r*c,min,max)
  dim(x)=c(r,c)
  return(x)
}

AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}

X_sim=function(r,c,group,min,max,tau){
  x=rnorm(r*c)
  dim(x)=c(r,c)
  covm=rep(0,c*c)
  dim(covm)=c(c,c)
  ng=length(group)
  for (i in 1:(ng-1)){
    ind=(group[i]+1):group[i+1]
    covm[ind,ind]=AR1(tau,group[i+1]-group[i])
    if(i!=(ng-1)){
      for (j in (i+1):(ng-1)) {
        ind2=(group[j]+1):group[j+1]
        cov_inter=Cov(group[i+1]-group[i],group[j+1]-group[j],min,max)
        covm[ind,ind2]=cov_inter
        covm[ind2,ind]=t(cov_inter)
      }
    }
    
  }
  return(X_gen(covm,x))
  
}

Simulation=function(n,p,Beta_true,snr){
  setwd("Desktop/research/Kevin He/PRS/Rewrite_Lassosum/BRIGHT/")
  library(Rcpp)
  library(mvtnorm)
  source("R/Std_preproc.R")
  dyn.load("C/my.so")
  sourceCpp("src/gdfit_Gaussian.cpp")
  
  n=10000
  p=100
  Beta_true=rep(0,p)
  Beta_true[1:10]=1:10
  snr=0.9
  
  X=rmvnorm(n,mean = rep(0,p),sigma = diag(rep(1,p)))
  Xbeta=X%*%Beta_true
  var_epsilon=var(Xbeta)*(1-snr)/snr
  Y=Xbeta+rnorm(length(Xbeta),0,sqrt(var_epsilon))
  
  
  XG <- newXG(X, 1:ncol(X), rep(1,ncol(X)), 1, FALSE)
  YG <- newY(Y,"gaussian")
  n <- nrow(XG$X)
  p <- ncol(XG$X)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  
  XtY=t(XG$X)%*%YG/n
  
  lambda.max=MaxLambda(XtY, tilde_beta=t(t(Beta_true)), XG$X, t(t(K1)), m=t(t(rep(1,p))), K0, tau=0,
            eta=0, alpha=1, eps=0.0001,max_iter=1000)
  lambda.min=0.001
  nlambda=100
  lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  
  rst=gdfit_gaussian(XtY, tilde_beta=t(t(Beta_true)), XG$X, t(t(lambda)), t(t(K1)), m=t(t(rep(1,p))), K0, penalty=1, tau=0,eta=0, alpha=1, gamma=0,
                 eps=0.0001,max_iter=1000, dfmax=p, gmax=p, user=T)
  
}
