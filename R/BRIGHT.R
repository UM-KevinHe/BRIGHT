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
Preprocess=function(chr, BP, A1, A2, p, n, sign,ref="EUR",group=1:p,m=NA){
  
  # read in reference panel
  if(ref%in%c("EUR","AFR","EAS","AMR","SAS")){
    KG=read.table(paste("inst/1000G/",ref,"_1000G_2504.bim",sep=""))
    gene_bed = BEDMatrix(paste("inst/1000G/",ref,"_1000G_2504.bed",sep=""))
    gene_bed=as.matrix(gene_bed)
  }
  else{
    KG=read.table(ref)
    gene_bed = BEDMatrix(ref)
    gene_bed=as.matrix(gene_bed)
  }
  
  #QC
  p[p==0]=10^-100
  XtY=p2cor(p,n,sign)
  pos_na=is.na(XtY)
  XtY=XtY[!pos_na]
  chr=chr[!pos_na]
  BP=BP[!pos_na]
  A1=A1[!pos_na]
  A2=A2[!pos_na]
  group=group[!pos_na]
  nm=paste(chr,BP,A2,A1,"b37",sep="_")
  nm_flp=paste(chr,BP,A1,A2,"b37",sep="_")
  
  pos_KG=as.logical(KG[,2]%in%nm+KG[,2]%in%nm_flp)
  KG=KG[pos_KG,]
  gene_bed=gene_bed[,pos_KG]
  pos_flp=KG[,2]%in%nm_flp
  gene_bed[,pos_flp]=2-gene_bed[,pos_flp]
  KG[pos_flp,2]=paste(KG[,1],KG[,4],KG[,5],KG[,6],"b37",sep="_")
  colnames(gene_bed)=KG[,2]
  
  pos_GWAS=as.logical(nm%in%KG[,2])
  XtY=XtY[pos_GWAS]
  chr=chr[pos_GWAS]
  BP=BP[pos_GWAS]
  A1=A1[pos_GWAS]
  A2=A2[pos_GWAS]
  group=group[pos_GWAS]
  rst=list()
  
  XG <- newXG(gene_bed, group, m, 1, FALSE)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  XtY <- t(t(cor[XG$nz]))
  
  if(is.na(sum(m))){
    m=sqrt(table(group[group!=0]))
  }
  
  LDB=read.table(paste("inst/data/Berisa.",ref,".hg19.bed",sep=""),header=T)
  LDB[,1]=as.numeric(sapply(strsplit(LDB[,1],"chr"), function(x){x[2]}))
  LDB=as.matrix(LDB)
  blk=t(t(c(LD(t(t(chr)),t(t(BP)),LDB),length(BP))))
  
  rst["XtY"]=XtY
  rst["XG"]=XG
  rst["m"]=m
  rst["blk"]=blk
  rst["K"]=K
  rst["K0"]=K0
  rst["K1"]=K1
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
BRIGHTs=function(XtY,Beta_prior,X,lambda,K1,m,blk,K0,penalty=1,tau=0,eta,alpha=1, gamma=0,eps=0.0000001,max_iter=1000000, dfmax=5000, gmax=5000, user=T){
  
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
