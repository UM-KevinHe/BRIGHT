setwd("Desktop/research/Kevin He/PRS/Rewrite_Lassosum/BRIGHT/")
library(Rcpp)
library(mvtnorm)
source("R/Std_preproc.R")
dyn.load("C/my.so")
sourceCpp("src/gdfit_Gaussian.cpp")

n=10000
p=100
Beta_true=rep(0,p)
Beta_true[1:20]=c(1:10,-(1:10))
snr=0.9

cor=0
sigma=diag(rep(1,p))
for (i in 2:p) {
  sigma[i,i-1]=cor
  sigma[i-1,i]=cor
}

X=rmvnorm(n,mean = rep(0,p),sigma = sigma)
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
eta=0
tilde_beta=Beta_true

lambda.max=MaxLambda(XtY, tilde_beta=t(t(tilde_beta)), XG$X, t(t(K1)), m=t(t(rep(1,p))), K0, tau=0,
                     eta=eta, alpha=1, eps=0.0001,max_iter=1000)
lambda.min=0.001
nlambda=100
lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))

lamb=lambda[52]

AIC=c()
TR=c()
for (k in 1:100){
  print(k)
  n=10000
  p=100
  Beta_true=rep(0,p)
  Beta_true[1:20]=c(1:10,-(1:10))
  snr=0.9
  
  cor=0
  sigma=diag(rep(1,p))
  for (i in 2:p) {
    sigma[i,i-1]=cor
    sigma[i-1,i]=cor
  }
  
  X=rmvnorm(n,mean = rep(0,p),sigma = sigma)
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
  eta=0
  tilde_beta=Beta_true
  
  X=rmvnorm(n*10,mean = rep(0,p),sigma = sigma)
  XG <- newXG(X, 1:ncol(X), rep(1,ncol(X)), 1, FALSE)
  
  
  rst=gdfit_gaussian(XtY, tilde_beta=t(t(tilde_beta)), XG$X, t(t(lamb)), t(t(K1)), m=t(t(rep(1,p))), K0, penalty=1, tau=0,eta=eta, alpha=1, gamma=0,
                     eps=0.0001,max_iter=100000, dfmax=p, gmax=p, user=T)
  
  OBJ=diag(as.matrix((1+eta)*t(rst$Beta)%*%rst$Sig%*%rst$Beta/2))-t(rst$Beta)%*%(XtY+eta*rst$Sig%*%tilde_beta)
  AIC=c(AIC,as.vector(OBJ+colSums(rst$Beta!=0)))
  ind=which.min(AIC)
  
  TR=c(TR,(t(rst$Beta-Beta_true)%*%sigma%*%(rst$Beta-Beta_true)+var_epsilon-Beta_true%*%sigma%*%Beta_true)/2)
}

# calculate MSE
TR_ind=c()
TR_sum=c()
AIC_ind=c()
AIC_sum=c()
df_ind=c()
df_sum=c()

for (k in 1:100){
  print(k)
  n=10000
  p=100
  Beta_true=rep(0,p)
  Beta_true[1:20]=c(1:10,-(1:10))
  snr=0.9
  
  cor=0
  sigma=diag(rep(1,p))
  for (i in 2:p) {
    sigma[i,i-1]=cor
    sigma[i-1,i]=cor
  }
  
  X=rmvnorm(n,mean = rep(0,p),sigma = sigma)
  Xbeta=X%*%Beta_true
  #var_epsilon=var(Xbeta)*(1-snr)/snr
  Y=Xbeta+rnorm(length(Xbeta),0,sqrt(var_epsilon))
  
  XG <- newXG(X, 1:ncol(X), rep(1,ncol(X)), 1, FALSE)
  YG <- newY(Y,"gaussian")
  n <- nrow(XG$X)
  p <- ncol(XG$X)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  
  XtY=t(XG$X)%*%YG/n
  eta=0
  tilde_beta=Beta_true
  
  rst=gdfit_gaussian(XtY, tilde_beta=t(t(tilde_beta)), XG$X, t(t(lambda)), t(t(K1)), m=t(t(rep(1,p))), K0, penalty=1, tau=0,eta=eta, alpha=1, gamma=0,
                         eps=0.0001,max_iter=100000, dfmax=p, gmax=p, user=T)
  
  OBJ=diag(as.matrix((1+eta)*t(rst$Beta)%*%rst$Sig%*%rst$Beta/2))-t(rst$Beta)%*%(XtY+eta*rst$Sig%*%tilde_beta)
  AIC_ind=rbind(AIC_ind,as.vector(OBJ+colSums(rst$Beta!=0)))
  df_ind=rbind(df_ind,as.vector(colSums(rst$Beta!=0)))
  tmp=c()
  for (a in 1:100) {
    tmp=c(tmp,(t(rst$Beta[,a]-Beta_true)%*%sigma%*%(rst$Beta[,a]-Beta_true)+var_epsilon-Beta_true%*%sigma%*%Beta_true)/2)
  }
  TR_ind=rbind(TR_ind,tmp)
  
  X=rmvnorm(n,mean = rep(0,p),sigma = sigma)
  XG <- newXG(X, 1:ncol(X), rep(1,ncol(X)), 1, FALSE)
  
  
  rst=gdfit_gaussian(XtY, tilde_beta=t(t(tilde_beta)), XG$X, t(t(lambda)), t(t(K1)), m=t(t(rep(1,p))), K0, penalty=1, tau=0,eta=eta, alpha=1, gamma=0,
                         eps=0.0001,max_iter=100000, dfmax=p, gmax=p, user=T)
  
  OBJ=diag(as.matrix((1+eta)*t(rst$Beta)%*%rst$Sig%*%rst$Beta/2))-t(rst$Beta)%*%(XtY+eta*rst$Sig%*%tilde_beta)
  AIC_sum=rbind(AIC_sum,as.vector(OBJ+colSums(rst$Beta!=0)))
  df_sum=rbind(df_sum,as.vector(colSums(rst$Beta!=0)))
  tmp=c()
  for (a in 1:100) {
    tmp=c(tmp,(t(rst$Beta[,a]-Beta_true)%*%sigma%*%(rst$Beta[,a]-Beta_true)+var_epsilon-Beta_true%*%sigma%*%Beta_true)/2)
  }
  TR_sum=rbind(TR_sum,tmp)
}

pdf("Ind_AIC_TR.pdf",width = 5,height = 5)
par(mar = c(4, 4, 0.2, 0.2))
plot(1,type = "n",ylim = c(-400,40),xlim = c(min(log(lambda)),max(log(lambda))),xlab = expression(log(lambda)),ylab = "AIC or TR",frame.plot = FALSE,cex.lab=1,cex.axis=1)
lines(log(lambda),colMeans(AIC_ind),ty="l",col="blue",lty=1,lwd=3)
lines(log(lambda),colMeans(TR_ind),ty="l",col="red",lty=,lwd=3)
legend("bottomright", legend=c("TR","AIC"),col=c("red", "blue"),lty=c(1,1),lwd=c(3,3))
dev.off()

pdf("sum_AIC_TR.pdf",width = 5,height = 5)
par(mar = c(4, 4, 0.2, 0.2))
plot(1,type = "n",ylim = c(-400,40),xlim = c(min(log(lambda)),max(log(lambda))),xlab = expression(log(lambda)),ylab = "AIC or TR",frame.plot = FALSE,cex.lab=1,cex.axis=1)
lines(log(lambda),colMeans(AIC_sum),ty="l",col="blue",lty=1,lwd=3)
lines(log(lambda),colMeans(TR_sum),ty="l",col="red",lty=,lwd=3)
legend("bottomright", legend=c("TR","AIC"),col=c("red", "blue"),lty=c(1,1),lwd=c(3,3))
dev.off()

pdf("Ind_AIC_TR_split.pdf",width = 5,height = 5)
par(mar = c(4, 4, 0.2, 0.2))
plot(1,type = "n",ylim = c(-400,40),xlim = c(min(log(lambda)),max(log(lambda))),xlab = expression(log(lambda)),ylab = "AIC or TR",frame.plot = FALSE,cex.lab=1,cex.axis=1)
for (k in 1:100) {
  lines(log(lambda),AIC_ind[k,],ty="l",col=alpha("blue",0.1),lty=1,lwd=1)
  lines(log(lambda),TR_ind[k,],ty="l",col=alpha("red",0.1),lty=1,lwd=1)
}
legend("bottomright", legend=c("TR","AIC"),col=c("red", "blue"),lty=c(1,1),lwd=c(3,3))
dev.off()

pdf("sum_AIC_TR_split.pdf",width = 5,height = 5)
par(mar = c(4, 4, 0.2, 0.2))
plot(1,type = "n",ylim = c(-400,40),xlim = c(min(log(lambda)),max(log(lambda))),xlab = expression(log(lambda)),ylab = "AIC or TR",frame.plot = FALSE,cex.lab=1,cex.axis=1)
for (k in 1:100) {
  lines(log(lambda),AIC_sum[k,],ty="l",col=alpha("blue",0.1),lty=1,lwd=1)
  lines(log(lambda),TR_sum[k,],ty="l",col=alpha("red",0.1),lty=1,lwd=1)
}
legend("bottomright", legend=c("TR","AIC"),col=c("red", "blue"),lty=c(1,1),lwd=c(3,3))
dev.off()

sz_AIC_ind=c()
sz_AIC_sum=c()
sz_TR_ind=c()
sz_TR_sum=c()
for(k in 1:100){
  sz_AIC_sum=c(sz_AIC_sum,df_sum[k,][which.min(AIC_sum[k,])])
  sz_AIC_ind=c(sz_AIC_ind,df_ind[k,][which.min(AIC_ind[k,])])
  sz_TR_sum=c(sz_TR_sum,df_sum[k,][which.min(TR_sum[k,])])
  sz_TR_ind=c(sz_TR_ind,df_ind[k,][which.min(TR_ind[k,])])
}

data=c(sz_AIC_sum,sz_AIC_ind,sz_TR_sum,sz_TR_ind)
data=data.frame(df=data,method=rep(c("AIC_sum","AIC_ind","TR_sum","TR_ind"),each=100))

pdf("df.pdf",width = 5,height = 5)
plot=ggplot(data, aes(x=method, y=df)) +
  geom_boxplot()+
  scale_fill_manual(values=c(2:7))+
  ylim(0,0.4)+
  theme_bw()+
  theme(axis.text=element_text(size=15,angle = 45,vjust = 0.5, hjust=1),
        axis.title=element_text(size=20),
        axis.title.x=element_blank(),
        legend.position = "none",)+
  ylab("R2")
print(plot)
dev.off()

pdf("solu_LasKL.pdf",width = 5,height = 5)
par(mar = c(4, 4, 0, 0))
plot(1,type = "n",ylim = c(-1,1),xlim = c(-log(0.1),-log(0.0013)),xlab = expression(-log(lambda)),ylab = "Estimated effect sizes",frame.plot = FALSE,cex.lab=1,cex.axis=1)
for (i in 1:length(SNP)) {
  if (!SNP[i]%in%Tr_eff_name ){
    lines(-seq(log(0.0013), log(0.1), length.out=40),rst_freq[i,1:40],ty="l",col=alpha("Grey",0.1),lty=1)}
}
for (i in 1:length(SNP)) {
  if (SNP[i]%in%Tr_eff_name ){
    lines(-seq(log(0.0013), log(0.1), length.out=40),rst_freq[i,1:40],ty="l",col="red")}
}
dev.off()


