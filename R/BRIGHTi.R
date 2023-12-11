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
BRIGHTi=function(data,m=NA,group=NA,Beta_prior,penalty=1,family="gaussian",eta=c(0,exp(seq(log(0.1),log(10),length.out=20))),lambda=NA,lambda.min=0.001,nlambda=100,tau=1/3,alpha=1, gamma=ifelse(penalty == 1, 4, 3),eps=0.01,max_iter=1000000, dfmax=5000, gmax=5000,user=T){
  KG=data$KG
  gene_bed = BEDMatrix(data$geno_dir)
  gene_bed=as.matrix(gene_bed)
  Y=data$Tpheno[,3]
  
  p=ncol(gene_bed)
  
  if(is.na(group)[1]){
    group=1:p
  }
  if(is.na(sum(m))){
    m=sqrt(table(group[group!=0]))
  }
  
  if(is.null(data$Pss$Coef)){
    warning("dat$Pss$Coef not found avoiding the usage of prior data")
    Beta_prior=rep(0,p)
    eta=0
  }else{
    Beta_prior=data$Pss$Coef
  }
  
  #yy <- newY(Y, family)
  yy <- (Y-mean(Y))/sd(Y)
  XG <- newXG(gene_bed, group, m, 1, FALSE)
  K <- as.integer(table(XG$g))
  K0 <- as.integer(if (min(XG$g)==0) K[1] else 0)
  K1 <- as.integer(if (min(XG$g)==0) cumsum(K) else c(0, cumsum(K)))
  yy=t(t(yy))
  if(is.na(lambda)){
    #lambda <- exp(seq(from=log(max(as.vector(cov(yy,XG$X)))),to=log(lambda.min),length.out=nlambda))
    Lmd=MaxLambdai(yy, tilde_beta=t(t(Beta_prior)), XG$X, t(t(K1)), m=t(t(rep(1,p))), K0, tau=tau, eta=0, alpha=alpha, eps=eps,max_iter=max_iter)
    lambda.max=Lmd$lambda.max
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))  
  }
  
  
  if(eta[1]!=0){
    eta_vec=c(0,eta)
  }else{
    eta_vec=eta
  }
  out=list()
  AIC=c()
  out[["lambda"]]=lambda
  out[["eta_vec"]]=eta_vec
  out[["Prior"]]=Beta_prior
  
  for(eta in eta_vec){
    rst_BRIGHTi=gdfit_gaussiani(yy, XG$X, t(t(Beta_prior)), t(t(lambda)), t(t(K1)), m=t(t(m)), K0, penalty, tau,eta, alpha, gamma,eps,max_iter, dfmax, gmax, user)
    AIC=cbind(AIC,2*(rst_BRIGHTi$dev2*nrow(XG$X)+colSums(as.matrix(rst_BRIGHTi$Beta)!=0)/(1+eta)))
    out[[as.character(eta)]]=rst_BRIGHTi
  }
  dinom=dim(AIC)[1]
  ind=which.min(AIC)
  rw=ind%%dinom
  cl=ind%/%dinom+1
  if(rw==0){
    rw=dinom
    cl=cl-1
  }
  Best_eta_AIC=eta_vec[cl]
  Best_lambda_AIC=out[["lambda"]][rw]
  Beta_BRIGHT_AIC=out[[as.character(Best_eta_AIC)]]$Beta[,rw]
  Beta_local_AIC=out[[as.character(0)]]$Beta[,which.min(AIC[,1])]
  
  if(length(eta_vec)==1){
    writeLines(strwrap(paste("Best lambda based on AIC is",round(Best_lambda_AIC,digits = 3))))
  }else{
    writeLines(strwrap(paste("Best eta based on AIC is",round(Best_eta_AIC,digits = 3))))
    writeLines(strwrap(paste("Best lambda based on AIC is",round(Best_lambda_AIC,digits = 3))))
  }
  
  out[["AIC"]]=AIC
  out[["Best_eta_AIC"]]=Best_eta_AIC
  out[["Best_lambda_AIC"]]=Best_lambda_AIC
  out[["Beta_BRIGHT_AIC"]]=Beta_BRIGHT_AIC
  out[["Beta_local_AIC"]]=Beta_local_AIC
  
  return(out)
}
