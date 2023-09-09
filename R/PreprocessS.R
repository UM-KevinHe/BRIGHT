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
PreprocessS=function(Tss, Tind, Tref, Pss=NA, Pind=NA, Pref=NA, TLDblocks="EUR", PLDblocks=NA,TN=NA,PN=NA){
  
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
  Tss=QC(Tss,Tind,KG,lab="target")
  
  if(sum(is.na(Pind))){
    rst=list()
    writeLines(strwrap("No prior summary statistics detected"))
    rst[["TLDblocks"]]=TLDblocks
    rst[["PLDblocks"]]=PLDblocks
    rst[["LD_dir"]]=LD_dir
    rst[["Tss"]]=Tss
    rst[["Tref"]]=Tref
    rst[["TN"]]=TN
    rst[["KG"]]=KG
    return(rst)
  }
  
  if(is.data.frame(Pss)){
    Pss=QC(Pss,Pind,KG,lab = "prior")
    Pind=colnames(Pss)[6]
    if(Pind!="Coef"){
      writeLines(strwrap('"Coef" not found for prior data fitting BRIGHTs with eta=0 for prior data and using AIC to select hyper-parameters'))
      Lsum=LASSOsum(Pss,Pind,PLDblocks,KG,PN,Pref,lab="Prior LASSOsum")
      Pss[,6]=Lsum$model$`0`$Beta[,which.min(Lsum$AIC)]
      colnames(Pss)[6]="Coef"
      #LASAIC=min(Lsum$AIC)
    }
  }else if(is.list(Pss)){
    for (i in 1:nrow(Pind)){
      Pss[[i]]=QC(Pss[[i]],Pind[i],KG,lab = paste("Prior",i))
    }
  }else {
    stop("Pss must be either a data frame of one prior data or a list containing multiple prior data")
  }
  rst=list()
  
  rst[["TLDblocks"]]=TLDblocks
  rst[["PLDblocks"]]=PLDblocks
  rst[["LD_dir"]]=LD_dir
  rst[["Tss"]]=Tss
  rst[["Pss"]]=Pss
  rst[["Tref"]]=Tref
  rst[["Pref"]]=Pref
  rst[["TN"]]=TN
  rst[["PN"]]=PN
  rst[["KG"]]=KG
  return(rst)
}
