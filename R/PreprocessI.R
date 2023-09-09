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

PreprocessI=function(Tgeno, Tpheno, Pss=NA, Pind=NA, Pref=NA, PLDblocks=NA,PN=NA){
  
  # read in genotype
  KG=read.table(paste(Tgeno,".bim",sep=""))
  KG[,2]=paste(KG[,1],KG[,4],KG[,5],KG[,6],sep="_")
  geno_dir=paste(Tgeno,".bed",sep="")
  
  writeLines(strwrap(paste(length(Pind)," source of prior information identified",sep="")))
  
  if(sum(is.na(Pind))){
    rst=list()
    writeLines(strwrap("No prior summary statistics detected"))
    rst[["Tgeno"]]=Tgeno
    rst[["Tpheno"]]=Tpheno
    rst[["KG"]]=KG
    return(rst)
  }
  
  if(is.data.frame(Pss)){
    Pss=QC(Pss,Pind,KG,lab = "Prior")
    Pind=colnames(Pss)[6]
    Lsum=LASSOsum(Pss,Pind,PLDblocks,KG,PN,Pref)
    Pss[,6]=Lsum$model$`0`$Beta[,which.min(Lsum$AIC)]
    colnames(Pss)[6]="Coef"
  }else if(is.list(Pss)){
    for (i in 1:nrow(Pind)){
      Pss[[i]]=QC(Pss[[i]],Pind[i],KG,lab = paste("Prior",i))
    }
  }else {
    stop("Pss must be either a data frame of one prior data or a list containing multiple prior data")
  }
  
  rst=list()
  rst[["Tgeno"]]=Tgeno
  rst[["Tpheno"]]=Tpheno
  rst[["geno_dir"]]=geno_dir
  rst[["PLDblocks"]]=PLDblocks
  rst[["Pss"]]=Pss
  rst[["Pref"]]=Pref
  rst[["PN"]]=PN
  rst[["KG"]]=KG
  return(rst)
}
