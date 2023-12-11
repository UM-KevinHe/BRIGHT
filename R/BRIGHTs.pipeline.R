#' Preprocess for summary statistics (BRIGHTs)
#'
#' Preprocess the GWAS and the reference. Quality control is performed first for both target and prior (please see QC section). Then prior "GWAS" or "Correlations" is converted to "Coefficients" through BRIGHTs with eta=0; the hyperparameters are tuned through BRCp.
#'
#' @param Tss a data frame of the target summary statistics, including "GWAS" or "Corr" corresponding to the GWAS or correlation summary statistics, the input summary type must match with Tind below.
#' @param Tind a character scaler of the input target summary type, it must be either "GWAS" or "Corr" corresponding to the GWAS or correlation summary statistics. It must match with the summary statistics in Tss.
#' @param Tref a character scaler indicating the LD reference. It takes values as either a directory to user provided Plink 1 binary files (.bim, .bed, .fam) or the 1000 genome project reference provided by the package: "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19). For the user provided reference, the SNPs must match with those provided in Pref with the same order.
#' @param Pss a data frame of the prior summary statistics, including "GWAS", "Corr" or "Coef" corresponding to the GWAS, correlation or coefficient summary statistics, the input summary type must match with Pind below.
#' @param Pind a character scaler of the input prior summary type, it must be either "GWAS", "Corr" or "Coef" corresponding to the GWAS, correlation or coefficient summary statistics. It must match with the summary statistics in Pss.
#' @param Pref a character scaler indicating the LD reference. It takes values as either a directory to user provided Plink 1 binary files (.bim, .bed, .fam) or the 1000 genome project reference provided by the package: "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19). For the user provided reference, the SNPs must match with those provided in Tref with the same order.
#' @param TLDblocks a character scaler indicating the LD block structure to be used from Berisa et al for the target population. It takes values from "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19). The genetic build must match with those provided in Tref.
#' @param PLDblocks a character scaler indicating the LD block structure to be used from Berisa et al for the prior population. It takes values from "EUR.hg19", "EUR.hg38", "AFR.hg19", "AFR.hg38", "EAS.hg19", "EAS.hg38", "AMR.hg19", "AMR.hg38", "SAS.hg19" or "SAS.hg38". The abbreviations before dot are populations (more detials can be found here: https://useast.ensembl.org/Help/Faq?id=532); The abbreviations after dot are the genetic build (more details can be found: https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19). The genetic build must match with those provided in Pref.
#' @param TN a integer scaler indicating the sample sizes of the target summary statistics input in Tss.
#' @param PN a integer scaler indicating the sample sizes of the target summary statistics input in Tss.
#'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
#'@param group a vector describing the grouping of the coefficients. For greatest efficiency and least ambiguity (see details), it is best if group is a factor or vector of consecutive integers, although unordered groups and character vectors are also allowed. If there are coefficients to be included in the model without being penalized, assign them to group 0 (or "0").
#'@param penalty a integer scalar for chosing the penalties; 1 corresponds to elastic net (with alpha=0 to be LASSO); 2 corresponds to MCP; 3 corresponds to SCAD.
#'@param tau tuning parameter for the group exponential lasso; defaults to 1/3.
#'@param eta a numeric scalar of the weights for incorporating prior information.
#'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
#'@param nlambda a integer scaler of the number of lambdas to be generated automatically (must be provided when lambda is not specified).
#'@param lambda.min a numeric scaler of the minimum lambda to for automatic generation (must be provided when lambda is not specified). 
#' @param alpha \code{BRIGHT} allows for both a group penalty and an L2 (ridge)
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
#' @return A list of preprocessed data that can be directly input to BRIGHTs().
BRIGHTs.pipeline = function(Tss, Tind, Tref, Pss = NA, Pind = NA, Pref = NA, TLDblocks = "EUR.hg19", PLDblocks = NA, TN = NA,
                            PN = NA, m = NA, group = NA, penalty = 1, tau = 1/3, eta = c(0,exp(seq(log(0.1),log(10), length.out = 20))),
                            lambda = NA, nlambda = 100, lambda.min = 0.001, alpha = 1, gamma = 0, eps = 0.0000001, max_iter = 1000000,
                            dfmax = 5000 , gmax = 5000, user = T){
  dat=PreprocessS(Tss, Tind, Tref, Pss, Pind, Pref, TLDblocks = c("SAS.hg19"), PLDblocks = c("EUR.hg19"), TN = 100, PN = 11514)
  out <- BRIGHTs(data = dat, m, group, penalty, tau, eta,
                 lambda, nlambda, lambda.min, alpha, gamma, eps, max_iter,
                 dfmax, gmax, user)
  return(out)
}
