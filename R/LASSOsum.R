LASSOsum=function(ss, ind, LDblocks, KG, N ,ref, lab="LASSOsum"){
  if(names(ss)[6]!="Coef"){
    ss=QC(ss,ind,KG,lab = lab)
    #writeLines(strwrap("Coefficients not identified, starting fitting LASSOsum on prior data"))
    data=list()
    data[["Tss"]]=ss
    data[["Tind"]]=ind
    data[["TLDblocks"]]=LDblocks
    data[["Tref"]]=ref
    data[["LD_dir"]]=paste(ref,".bed",sep="")
    data[["TN"]]=N
    Prior_Lassosum=BRIGHTs(data,m=NA,group=NA,penalty=1,tau=0,eta=0,lambda=NA,nlambda=50,lambda.min=0.0001,alpha=1, gamma=0,eps=0.001,max_iter=1000000, dfmax=5000, gmax=5000, user=T)
    if(is.na(N)){
      writeLines(strwrap("sample size of prior data must be provided for the AIC calculation, otherwise please directly input Coefficients from joint models"))
    }else{
      AIC=2*(Prior_Lassosum[[as.character("0")]]$dev2*N+colSums(as.matrix(Prior_Lassosum[[as.character(0)]]$Beta)!=0))
    }
    rst=list()
    rst[["AIC"]]=AIC
    rst[["model"]]=Prior_Lassosum
    return(rst)
  }else{
    stop('"Coef" found but "IProd" or "GWAS required for LASSOsum fitting, aborting"')
  }
}
