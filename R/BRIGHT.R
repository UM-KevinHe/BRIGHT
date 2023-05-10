#' GWAS to correlation
#'
#' Recover marginal correlation information from GWAS summary statistics.
#'
#' @param p a numeric vector of the p values from marginal SLR.
#' @param n a numeric vector of the sample sizes.
#' @param sign a numeric vector of the signs of coefficient estimates.
#' @return Recovered XtY/n from the GWAS summary.
#' @examples 
#' p=c(0.01,0.015)
#' n=c(100,1000)
#' sign=c(-1,1)
#' p2cor(p,n,sign)
p2cor <- function(p, n, sign=rep(1, length(p))) {
  
  stopifnot(length(n)==1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))
  
  t <- sign(sign) * qt(p/2, df=n-2, lower.tail=F)
  
  return(t / sqrt(n - 2 + t^2))
  
}

QC <- function(ss,ind,KG){
  if(ind=="GWAS"){
    ss$P[ss$P==0]=10^-100
    XtY=p2cor(ss$P,ss$N,ss$Sign)
    outnm="IProd"
  }else if (ind=="IProd"){
    XtY=ss$IProd
    outnm="IProd"
  }else if (ind=="Coef"){
    XtY=ss$Coef
    outnm="Coef"
  }else{
    stop('The data type indicator (Tind or Pind) need to be one of "GWAS", "IProd", or "Coef"')
  }
  
  pos_na=is.na(XtY)
  XtY=XtY[!pos_na]
  ss=ss[!pos_na,]
  nm=paste(ss$CHR,ss$BP,ss$A1,ss$A2,sep="_")
  nm_flp=paste(ss$CHR,ss$BP,ss$A2,ss$A1,sep="_")

  New_ss=cbind(nm,nm_flp,ss[,1:4],XtY)
  colnames(New_ss)=c("nm","nm_flp","CHR","BP","A1","A2",outnm)
  
  pos_ss=as.logical(New_ss$nm%in%KG[,2]+New_ss$nm_flp%in%KG[,2])
  New_ss=New_ss[pos_ss,]
  pos_flp=New_ss$nm_flp%in%KG[,2]
  New_ss[pos_flp,1]=New_ss[pos_flp,2]
  temp_A1=New_ss[pos_flp,5]
  New_ss[pos_flp,5]=New_ss[pos_flp,6]
  New_ss[pos_flp,6]=temp_A1
  New_ss[pos_flp,7]=-New_ss[pos_flp,7]
  New_ss=New_ss[,-2]
  return(New_ss)
}


#' Preprocess
#'
#' Preprocess the GWAS and the 1000 genome reference.
#'
#' @param chr a numeric vector of the chromesomes of each SNP in GWAS.
#' @param BP a numeric vector of the position in base pair of each SNP in GWAS.
#' @param A1 a string vector of the effect allele of each SNP in GWAS.
#' @param A2 a string vector of the reference allele of each SNP in GWAS.
#' @param p a numeric vector of the p values from marginal t-test.
#' @param n a numeric vector of the sample sizes.
#' @param sign a numeric vector of the signs of coefficient estimates.
#' @param ref a string of the GWAS population; EUR, European ancestry; AFR, African ancestry; AMR, admixed American; EAS, East Asian; SAS, South Asian.
#' @param group a vector describing the grouping of the coefficients.  For
#' greatest efficiency and least ambiguity (see details), it is best if
#' \code{group} is a factor or vector of consecutive integers, although
#' unordered groups and character vectors are also allowed.  If there are
#' coefficients to be included in the model without being penalized, assign
#' them to group 0 (or \code{"0"}).
#' @param m a vector of values representing multiplicative
#' factors by which each group's penalty is to be multiplied.  Often, this is a
#' function (such as the square root) of the number of predictors in each
#' group.  The default is to use the square root of group size for the group
#' selection methods, and a vector of 1's (i.e., no adjustment for group size)
#' for bi-level selection.
#' @return A list of preprocessed XtY and reference genotype matrix.
PreprocessS=function(Tss, Tind, Tref, Pss, Pind, Pref, LDblocks="EUR"){
  
  # read in reference panel
  if(Tref%in%c("EUR","AFR","EAS","AMR","SAS")){
    writeLines(strwrap(paste("Using ",Tref,"LD from 1000 genome project",sep="")))
    KG=read.table(paste("inst/1000G/",Tref,"_1000G_2504.bim",sep=""))
    KG[,2]=paste(KG[,1],KG[,4],KG[,5],KG[,6],sep="_")
    LD_dir=paste("inst/1000G/",ref,"_1000G_2504.bed",sep="")
  }else{
    writeLines(strwrap("Using user specified LD"))
    KG=read.table(paste(Tref,".bim",sep=""))
    KG[,2]=paste(KG[,1],KG[,4],KG[,5],KG[,6],sep="_")
    LD_dir=paste(Tref,".bed",sep="")
  }
  
  writeLines(strwrap(paste(length(Pind)," source of prior information identified",sep="")))
  
  #QC
  Tss=QC(Tss,Tind,KG)
  writeLines(strwrap(paste(nrow(Tss)," SNPs passed filtering for target data",sep="")))
  
  if(is.data.frame(Pss)){
    Pss=QC(Pss,Pind,KG)
    writeLines(strwrap(paste(nrow(Pss)," SNPs passed filtering for prior data 1",sep="")))
    #Prior_Lassosum=BRIGHTs(Tss=Pss,Tref=Pref,LDblocks="EUR",Pss=NA,m=NA,group=NA,m=NA,penalty=1,tau=0,eta=0,lambda=NA,nlambda=100,lambda.min=0.001,alpha=1, gamma=0,eps=0.0001,max_iter=1000000, dfmax=5000, gmax=5000, user=T)
  }else if(is.list(Pss)){
    for (i in 1:nrow(Pind)){
      Pss[[i]]=QC(Pss[[i]],Pind[i],KG)
      writeLines(strwrap(paste(nrow(Pss[[i]])," SNPs passed filtering for prior data ",i,sep="")))
    }
  }else {
    stop("Pss must be either a data frame of one prior data or a list containing multiple prior data")
  }
  rst=list()
  
  rst[["LD_dir"]]=LD_dir
  rst[["Tss"]]=Tss
  rst[["Pss"]]=Pss
  rst[["Tref"]]=Tref
  rst[["Pref"]]=Pref
  rst[["KG"]]=KG
  return(rst)
}

#' BRIGHT estimation procedure for summary level data
#'
#' Fit BRIGHT estimation procedure for summary level data with grouped penalties over a grid of values of the regularization parameter lambda.
#'
#' For more information about the penalties and their properties, please
#' consult the references below, many of which contain discussion, case
#' studies, and simulation studies comparing the methods.  If you use
#' \code{BRIGHT} for an analysis, please cite the appropriate reference.
#' 
#' In keeping with the notation from the original MCP paper, the tuning
#' parameter of the MCP penalty is denoted 'gamma'.  Note, however, that in
#' Breheny and Huang (2009), \code{gamma} is denoted 'a'.
#' 
#' The objective function for \code{BRIGHTi} optimization is defined to be
#' \deqn{Q(\beta)=r\beta+\beta\Sigma\beta/2+(\beta'-\beta)\Sigma(\beta'-\beta)/2n+p(\beta),}  where the first two terms on the right hand side (RHS) are
#' a approximation of the OLS loss; the second term on the RHS is the Bregman-divergence; the third term on the RHS is the penalty.
#' 
#' The algorithms employed by \code{BRIGHT} are stable and generally converge
#' quite rapidly to values close to the solution.  However, especially when p
#' is large compared with n, \code{BRIGHT} may fail to converge at low values
#' of \code{lambda}, where models are nonidentifiable or nearly singular.
#' Often, this is not the region of the coefficient path that is most
#' interesting.  The default behavior warning the user when convergence
#' criteria are not met may be distracting in these cases, and can be modified
#' with \code{warn} (convergence can always be checked later by inspecting the
#' value of \code{iter}).
#' 
#' If models are not converging, increasing \code{max.iter} may not be the most
#' efficient way to correct this problem.  Consider increasing \code{n.lambda}
#' or \code{lambda.min} in addition to increasing \code{max.iter}.
#' 
#' Although \code{BRIGHT} allows groups to be unordered and given arbitary
#' names, it is recommended that you specify groups as consecutive integers.
#' The first reason is efficiency: if groups are out of order, \code{X} must be
#' reordered prior to fitting, then this process reversed to return
#' coefficients according to the original order of \code{X}.  This is
#' inefficient if \code{X} is very large.  The second reason is ambiguity with
#' respect to other arguments such as \code{group.multiplier}.  With
#' consecutive integers, \code{group=3} unambiguously denotes the third element
#' of \code{group.multiplier}.
#' 
#' \code{BRIGHT} requires groups to be non-overlapping.
#'
#'@param XtY a numeric verticle vector of the marginal correlation between genotype and outcome, can be recovered from GWAS summary statistics by function p2cor().
#'@param Beta_prior a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
#'@param X a numeric matrix of reference genotypes, default is from 1000 genome project.
#'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
#'@param K1 a numeric verticle vector of the grouping of SNPs for group penalties; SNPs within K1[i,i+1] are grouped together.
#'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
#'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
#'@param K0 a integer scalar for the number of SNPs/covariates that is not penalized.
#'@param penalty a integer scalar for chosing the penalties; 1 corresponds to LASSO, 2 corresponds to MCP, 3 corresponds to SCAD.
#' @param tau tuning parameter for the group exponential lasso; defaults to
#' 1/3.
#'@param eta a numeric scalar of the weights for incorporating prior information.
#' @param alpha \code{grpreg} allows for both a group penalty and an L2 (ridge)
#' penalty; \code{alpha} controls the proportional weight of the regularization
#' parameters of these two penalties.  The group penalties' regularization
#' parameter is \code{lambda*alpha}, while the regularization parameter of the
#' ridge penalty is \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.
#' @param gamma tuning parameter of the group or composite MCP/SCAD penalty
#' (see details).  Default is 3 for MCP and 4 for SCAD.
#' @param eps convergence threshhold.  The algorithm iterates until the RMSD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{1e-4}.  See details.
#' @param max.iter maximum number of iterations (total across entire path).
#' Default is 10000.  See details.
#' @param dfmax limit on the number of parameters allowed to be nonzero.  If
#' this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gmax limit on the number of groups allowed to have nonzero elements.
#' If this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#'@param user a logical scalar for indicating whether lambda is user specified; if user=TRUE, then the iteration will start from the largest lambda, otherwise the iteration will start from the second largest lambda.
#'@return Write a list of results containing: beta, the coefficient estimate from BRIGHT; 
#'@return iter, the number of total iterations needed for the model to converge with each lambda; 
#'@return df total degree of freedom of the converged model with each lambda;
#'@return dev, the approximated deviance associated with each lambda.
BRIGHTs=function(Tss,Tref,LDblocks="EUR",Pss=NA,m=NA,group=NA,penalty=1,tau=0,eta=0,lambda=NA,nlambda=100,lambda.min=0.001,alpha=1, gamma=0,eps=0.0000001,max_iter=1000000, dfmax=5000, gmax=5000, user=T){
  p=nrow(Tss)
  group=1:p
  
  if(is.na(Pss)[1]){
    Beta_prior=rep(0,p)
    eta=0
  }else{
    Beta_prior=Pss$Coef
  }
  
  if(is.na(sum(m))){
    m=sqrt(table(group[group!=0]))
  }
  
  gene_bed = BEDMatrix(paste(Tref,".bed",sep=""))
  gene_bed=as.matrix(gene_bed)
  
  XG <- newXG(gene_bed, group, m, 1, FALSE)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  XtY <- t(t(Tss$IProd))
  
  LDB=read.table(paste("inst/data/Berisa.",LDblocks,".hg19.bed",sep=""),header=T)
  LDB[,1]=as.numeric(sapply(strsplit(as.vector(LDB[,1]),"chr"), function(x){x[2]}))
  LDB=as.matrix(LDB)
  blk=t(t(c(LD(t(t(Tss$CHR)),t(t(Tss$BP)),LDB),length(Tss$BP))))
  
  if(is.na(sum(lambda))){
    Lmd=MaxLambda(XtY, tilde_beta=t(t(Beta_prior)), XG$X, t(t(K1)), m=t(t(rep(1,p))), blk, K0, tau=tau,
                  eta=eta, alpha=alpha, eps=eps,max_iter=max_iter)
    lambda.max=Lmd$lambda.max
    lambda.min=0.001
    nlambda=100
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  }

  
  rst_BRIGHTs=gdfit_gaussian(XtY, t(t(Beta_prior)), XG$X, t(t(lambda)), t(t(K1)), m=t(t(m)), blk, K0, penalty, tau,eta, alpha, gamma,eps,max_iter, dfmax, gmax, user)

  return(rst_BRIGHTs)
}

#' BRIGHT estimation procedure for individual level data
#' 
#' Fit BRIGHT estimation procedure for individual level data with grouped penalties over a grid of values of the regularization parameter lambda.
#' 
#' For more information about the penalties and their properties, please
#' consult the references below, many of which contain discussion, case
#' studies, and simulation studies comparing the methods.  If you use
#' \code{BRIGHT} for an analysis, please cite the appropriate reference.
#' 
#' In keeping with the notation from the original MCP paper, the tuning
#' parameter of the MCP penalty is denoted 'gamma'.  Note, however, that in
#' Breheny and Huang (2009), \code{gamma} is denoted 'a'.
#' 
#' The objective function for \code{BRIGHTi} optimization is defined to be
#' \deqn{Q(\beta)=||y-X\beta||/2n+||X\beta'-X\beta||/2n+p(\beta),}  where the first term on the right hand side (RHS) is
#' the OLS loss; the second term on the RHS is the Bregman-divergence; the third term on the RHS is the penalty.
#' 
#' The algorithms employed by \code{BRIGHT} are stable and generally converge
#' quite rapidly to values close to the solution.  However, especially when p
#' is large compared with n, \code{BRIGHT} may fail to converge at low values
#' of \code{lambda}, where models are nonidentifiable or nearly singular.
#' Often, this is not the region of the coefficient path that is most
#' interesting.  The default behavior warning the user when convergence
#' criteria are not met may be distracting in these cases, and can be modified
#' with \code{warn} (convergence can always be checked later by inspecting the
#' value of \code{iter}).
#' 
#' If models are not converging, increasing \code{max.iter} may not be the most
#' efficient way to correct this problem.  Consider increasing \code{n.lambda}
#' or \code{lambda.min} in addition to increasing \code{max.iter}.
#' 
#' Although \code{BRIGHT} allows groups to be unordered and given arbitary
#' names, it is recommended that you specify groups as consecutive integers.
#' The first reason is efficiency: if groups are out of order, \code{X} must be
#' reordered prior to fitting, then this process reversed to return
#' coefficients according to the original order of \code{X}.  This is
#' inefficient if \code{X} is very large.  The second reason is ambiguity with
#' respect to other arguments such as \code{group.multiplier}.  With
#' consecutive integers, \code{group=3} unambiguously denotes the third element
#' of \code{group.multiplier}.
#' 
#' \code{BRIGHT} requires groups to be non-overlapping.
#' 
#' @param X_bed the genotype matrix directory stored as plink .bed format.  \code{BRIGHT} automatically load it and standardizes the data.
#' 
#' @param y the response vector.
#' @param group A vector describing the grouping of the coefficients.  For
#' greatest efficiency and least ambiguity (see details), it is best if
#' \code{group} is a factor or vector of consecutive integers, although
#' unordered groups and character vectors are also allowed.  If there are
#' coefficients to be included in the model without being penalized, assign
#' them to group 0 (or \code{"0"}).
#' @param Beta_prior a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
#' @param penalty the penalty to be applied to the model.  For group selection,
#' one of \code{grLasso}, \code{grMCP}, or \code{grSCAD}.
#' @param family either "gaussian" or "binomial", depending on the response.
#' @param eta a numeric scalar of the weights for incorporating prior information.
#' @param lambda a user supplied sequence of \code{lambda} values.  Typically,
#' this is left unspecified, and the function automatically computes a grid of
#' lambda values that ranges uniformly on the log scale over the relevant range
#' of lambda values.
#' @param lambda.min the smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}.  Default is .0001 if the number of observations is larger
#' than the number of covariates and .05 otherwise.
#' @param alpha \code{grpreg} allows for both a group penalty and an L2 (ridge)
#' penalty; \code{alpha} controls the proportional weight of the regularization
#' parameters of these two penalties.  The group penalties' regularization
#' parameter is \code{lambda*alpha}, while the regularization parameter of the
#' ridge penalty is \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.
#' @param eps convergence threshhold.  The algorithm iterates until the RMSD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{1e-4}.  See details.
#' @param max.iter maximum number of iterations (total across entire path).
#' Default is 10000.  See details.
#' @param dfmax limit on the number of parameters allowed to be nonzero.  If
#' this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gmax limit on the number of groups allowed to have nonzero elements.
#' If this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gamma tuning parameter of the group or composite MCP/SCAD penalty
#' (see details).  Default is 3 for MCP and 4 for SCAD.
#' @param tau tuning parameter for the group exponential lasso; defaults to
#' 1/3.
#' 
#' @return An object with S3 class \code{"grpreg"} containing:
#' \describe{
#' \item{beta}{The fitted matrix of coefficients.  The number of rows is equal
#' to the number of coefficients, and the number of columns is equal to
#' \code{nlambda}.}
#' \item{family}{Same as above.}
#' \item{group}{Same as above.}
#' \item{lambda}{The sequence of \code{lambda} values in the path.}
#' \item{alpha}{Same as above.}
#' \item{loss}{A vector containing either the residual sum of squares (`"gaussian"`) or negative log-likelihood (`"binomial"`) of the fitted model at each value of `lambda`.}
#' \item{n}{Number of observations.}
#' \item{penalty}{Same as above.}
#' \item{df}{A vector of length `nlambda` containing estimates of effective number of model parameters all the points along the regularization path.  For details on how this is calculated, see Breheny and Huang (2009).}
#' \item{iter}{A vector of length `nlambda` containing the number of iterations until convergence at each value of `lambda`.}
#' \item{group.multiplier}{A named vector containing the multiplicative constant applied to each group's penalty.}
#' }
#' 
#' @author Qinmengge Li
#' 
#' @seealso [BRIGHTs()] methods.
#' 
#' @references
#' \itemize{
#' \item Li Q, Patrick MT, Zhang H, Khunsriraksakul C, Stuart PE, Gudjonsson JE,
#'  Nair R, Elder JT, Liu DJ, Kang J, Tsoi LC. Bregman Divergence-Based Data 
#'  Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment.
#'   arXiv preprint arXiv:2210.06025. 2022 Oct 12.
#' 
#' \item Breheny P and Huang J. (2009) Penalized methods for bi-level variable
#' selection. *Statistics and its interface*, **2**: 369-380.
#' \doi{10.4310/sii.2009.v2.n3.a10}
#' 
#' \item Huang J, Breheny P, and Ma S. (2012). A selective review of group
#' selection in high dimensional models. *Statistical Science*, **27**: 481-499.
#' \doi{10.1214/12-sts392}
#' 
#' \item Breheny P and Huang J. (2015) Group descent algorithms for nonconvex
#' penalized linear and logistic regression models with grouped predictors.
#' *Statistics and Computing*, **25**: 173-187. \doi{10.1007/s11222-013-9424-2}
#' 
#' \item Breheny P. (2015) The group exponential lasso for bi-level variable
#' selection. *Biometrics*, **71**: 731-740. \doi{10.1111/biom.12300}
#' }
BRIGHTi=function(X_plink,Y,group=1:ncol(X),Beta_prior,penalty="grLasso",family="gaussian",eta,lambda,lambda.min,tau=1/3,alpha=1, gamma=ifelse(penalty == "grSCAD", 4, 3),eps=0.0000001,max_iter=1000000, dfmax=5000, gmax=5000){
  KG=read.table(paste(X_plink,".bim",sep=""))
  gene_bed = BEDMatrix(paste(X_plink,".bed",sep=""))
  gene_bed=as.matrix(gene_bed)
  
  rst_BRIGHTi=grpreg(X=gene_bed,y=Y,group=group,penalty = penalty,family = family,lambda = lambda,lambda.min = lambda.min,alpha=alpha,eps=eps,max.iter = max_iter,dfmax = dfmax,gmax = gmax,tau = tau,gamma = gamma)
  
  return(rst_BRIGHTi)
}

Valid.Ind <- function(out, Testpheno, Testgeno){
  gene_bed = BEDMatrix(paste(Testgeno,".bed",sep=""))
  gene_bed=as.matrix(gene_bed)
  XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
  XB_Exom=XG$X%*%out[["Prior"]]
  phe=Testpheno[,3]
  phe_std=(phe-mean(phe))/sd(phe)
  LAScor=cor(XB_Exom, phe_std)
  LASMSE=colMeans((XB_Exom-as.vector(phe_std))^2)
  eta_vec=out[["eta_vec"]]
  
  value=rep(0,6)
  pos=rep(0,6)
  eta_best=rep(0,3)
  ER=c()
  BRMSE_rst=c()
  BRCor_rst=c()
  Vali_MSE_rst=c()
  Vali_cor_rst=c()
  AIC_rst=c()
  Beta=list()
  Vali_MSE_best=Inf
  Vali_cor_best=-Inf
  AIC_best=Inf
  XB_MSE=c()
  XB_Cor=c()
  
  phe_test_std=(phe-mean(phe))/sd(phe)
  beta=out[[as.character(0)]]
  XB_Ind2=XG$X%*%beta
  Vali_MSE=colMeans((XB_Ind2-as.vector(phe_test_std))^2)
  Vali_cor=cor(XB_Ind2, phe_test_std)
  XB_local=XB_Ind2[,which.max(Vali_cor)]
  
  for(eta in eta_vec){
    beta=out[[as.character(eta)]]
    
    XB_Ind2=XG$X%*%beta
    
    Vali_MSE=colMeans((XB_Ind2-as.vector(phe_test_std))^2)
    Vali_cor=cor(XB_Ind2, phe_test_std)
    Vali_MSE_ind=which.min(Vali_MSE)
    Vali_cor_ind=which.max(Vali_cor)
    
    BRcor=cor(XB_Ind2, phe)
    phe_std=(phe-mean(phe))/sd(phe)
    BRMSE=colMeans((XB_Ind2-as.vector(phe_std))^2)
    
    if(min(Vali_MSE,na.rm=T)<Vali_MSE_best){
      Vali_MSE_best=min(Vali_MSE,na.rm=T)
      value[1]=BRcor[Vali_MSE_ind]
      value[2]=BRMSE[Vali_MSE_ind]
      eta_best[1]=eta
      XB_MSE=XB_Ind2[,which.min(Vali_MSE)]
    }
    if(max(Vali_cor,na.rm=T)>Vali_cor_best){
      Vali_cor_best=max(Vali_cor,na.rm=T)
      value[3]=BRcor[Vali_cor_ind]
      value[4]=BRMSE[Vali_cor_ind]
      eta_best[2]=eta
      XB_Cor=XB_Ind2[,which.max(Vali_cor)]
    }
    
    Vali_MSE_rst=cbind(Vali_MSE_rst,Vali_MSE)
    Vali_cor_rst=cbind(Vali_cor_rst,Vali_cor)
    
    BRMSE_rst=cbind(BRMSE_rst,BRMSE)
    BRCor_rst=cbind(BRCor_rst,BRcor)
    
  }
  rst=list()
  rst[["BRMSE_rst"]]=BRMSE_rst
  rst[["BRCor_rst"]]=BRCor_rst
  rst[["LAScor"]]=LAScor
  rst[["LASMSE"]]=LASMSE
  rst[["XB_Exom"]]=XB_Exom
  rst[["XB_local"]]=XB_local
  rst[["XB_MSE"]]=XB_MSE
  rst[["XB_Cor"]]=XB_Cor
  rst[["phe"]]=phe
  rst[["phe_std"]]=phe_std
  rst[["eta_vec"]]=eta_vec
  return(rst)
}

MSE_Cor.plot <- function(Val){
  eta_vec=Val[["eta_vec"]]
  dinom=dim(Val[["BRMSE_rst"]])[1]
  
  m=matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
  
  layout(mat = m,heights = c(0.1,0.9))
  
  par(mar = c(0, 0, 0, 0))
  plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
  legend(x="center", legend=c("BRIGHT","Local","Prior"),col=c("black", "black","black"),lty=c(1,NA,NA),lwd=c(3,NA,NA),pch = c(NA,15,16),cex = c(1.5,1.5,1.5),bty = "n",horiz = T)
  
  BRCor=Val[["BRMSE_rst"]]
  eta=eta_vec[-22]
  
  ind=which.min(BRCor)
  
  rw=ind%%dinom
  cl=ind%/%dinom+1
  ER=BRCor[,-22]
  
  ER_Exom=Val[["LASMSE"]]
  
  par(mar = c(5, 5, 0.5, 0.5))
  plot(1,type = "n",ylim = c(min(ER[rw,]),max(ER[rw,],ER_Exom)),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = "MSPE",cex.lab=1.5,cex.axis=1.5)
  lines(eta,ER[rw,],lwd=3)
  points(eta[1],ER[rw,1],pch=15,cex=2)
  points(eta[21],ER_Exom,pch=16,cex=2)
  
  BRCor=Val[["BRCor_rst"]]
  eta=eta_vec[-22]
  
  ind=which.max(BRCor)
  
  rw=ind%%dinom
  cl=ind%/%dinom+1
  ER=BRCor[,-22]^2
  
  ER_Exom=Val[["LAScor"]]^2
  
  par(mar = c(5, 5, 0.5, 0.5))
  plot(1,type = "n",ylim = c(min(ER[rw,],ER_Exom),max(ER[rw,],ER_Exom)),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = expression(R^2),cex.lab=1.5,cex.axis=1.5)
  lines(eta,ER[rw,],lwd=3)
  points(eta[1],ER[rw,1],pch=15,cex=2)
  points(eta[21],ER_Exom,pch=16,cex=2)
}

Density.plot <- function(Val,Pct){
  XB_local=Val[["XB_local"]]
  XB_prior=Val[["XB_Exom"]]
  XB_BRIGHT=Val[["XB_Cor"]]
  phe=Val[["phe"]]
  
  ind_local=XB_local>quantile(XB_local,probs = c(0,0.1,0.5,Pct,1))[4]
  ind_prior=XB_prior>quantile(XB_prior,probs = c(0,0.1,0.5,Pct,1))[4]
  ind_BRIGHT=XB_BRIGHT>quantile(XB_BRIGHT,probs = c(0,0.1,0.5,Pct,1))[4]
  
  high_local=phe[ind_local]
  low_local=phe[!ind_local]
  
  high_prior=phe[ind_prior]
  low_prior=phe[!ind_prior]
  
  high_BRIGHT=phe[ind_BRIGHT]
  low_BRIGHT=phe[!ind_BRIGHT]
  
  m=matrix(c(1,1,1,2,3,4),nrow = 2,ncol = 3,byrow = TRUE)
  
  layout(mat = m,heights = c(0.1,0.9))
  
  par(mar = c(0, 0, 0, 0))
  plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
  legend(x="center", legend=c(paste("Upper ",round(Pct*100)," PRS percentile",sep = ""),paste("Lower ",round(Pct*100)," PRS percentile",sep = "")),col=c("red", "black"),lty=c(2,1),lwd=c(3,3),cex = c(2,2),bty = "n",horiz = T)
  
  dist_local=density(high_local)$x[which.max(density(high_local)$y)]-density(low_local)$x[which.max(density(low_local)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.12),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("Local (",round(dist_local,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_local),lty=2,lwd=3,col="red")
  lines(density(low_local),lty=1,lwd=3)
  abline(v = density(high_local)$x[which.max(density(high_local)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_local)$x[which.max(density(low_local)$y)], col="black", lwd=3, lty=3)
  
  dist_prior=density(high_prior)$x[which.max(density(high_prior)$y)]-density(low_prior)$x[which.max(density(low_prior)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.12),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("Prior (",round(dist_prior,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_prior),lty=2,lwd=3,col="red")
  lines(density(low_prior),lty=1,lwd=3)
  abline(v = density(high_prior)$x[which.max(density(high_prior)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_prior)$x[which.max(density(low_prior)$y)], col="black", lwd=3, lty=3)
  
  dist_BRIGHT=density(high_BRIGHT)$x[which.max(density(high_BRIGHT)$y)]-density(low_BRIGHT)$x[which.max(density(low_BRIGHT)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.12),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("BRIGHT (",round(dist_BRIGHT,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_BRIGHT),lty=2,lwd=3,col="red")
  lines(density(low_BRIGHT),lty=1,lwd=3)
  abline(v = density(high_BRIGHT)$x[which.max(density(high_BRIGHT)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_BRIGHT)$x[which.max(density(low_BRIGHT)$y)], col="black", lwd=3, lty=3)
}

ROC.plot <- function(Val, Pct=0.5){
  XB_local=Val[["XB_local"]]
  XB_prior=Val[["XB_Exom"]]
  XB_BRIGHT=Val[["XB_Cor"]]
  phe=Val[["phe"]]
  
  phe_bi=phe>quantile(phe,probs = c(Pct))[1]
  
  suppressWarnings({
    ROC_local=pROC::roc(phe_bi~XB_local)
  })
  suppressWarnings({
    ROC_prior=pROC::roc(phe_bi~as.vector(XB_prior))
  })
  suppressWarnings({
    ROC_BRIGHT=pROC::roc(phe_bi~XB_BRIGHT)
  })
  
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,1),xlim = c(1,0),xlab = "Specificity",ylab = "Sensitivity",cex.lab=2,cex.axis=2,cex.main=2,frame.plot=F)
  box(bty="l")
  lines(ROC_BRIGHT,col="black",lty=1,lwd=3)
  lines(ROC_prior,col="red",lty=2,lwd=2)
  lines(ROC_local,col="blue",lty=3,lwd=2)
  legend("bottomright", legend=c(paste("BRIGHT (AUC:",round(as.numeric(ROC_BRIGHT$auc),digits = 3),")",sep = ""),paste("Prior (AUC:",round(as.numeric(ROC_prior$auc),digits = 3),")",sep = ""),paste("Local (AUC:",round(as.numeric(ROC_local$auc),digits = 3),")",sep = "")),col=c("black", "blue", "red"),lty=c(1,3,2), lwd = c(3,2,2), cex=1)
}