QC <- function(ss,ind,KG,lab){
  if(ind=="GWAS"){
    check=c("CHR","BP","A1","A2","Sign","N","P")
    if(sum(!check%in%colnames(ss))){
      stop(paste(check[!check%in%colnames(ss)]," column not found in ",lab," summary statistics",sep=""))
    }
    colnames(ss)
    ss$P[ss$P==0]=10^-100
    XtY=p2cor(ss$P,ss$N,ss$Sign)
    outnm="Corr"
  }else if (ind=="Corr"){
    check=c("CHR","BP","A1","A2","Corr")
    if(sum(!check%in%colnames(ss))){
      stop(paste(check[!check%in%colnames(ss)]," column not found in ",lab," summary statistics",sep=""))
    }
    XtY=ss$Corr
    outnm="Corr"
  }else if (ind=="Coef"){
    check=c("CHR","BP","A1","A2","Coef")
    if(sum(!check%in%colnames(ss))){
      stop(paste(check[!check%in%colnames(ss)]," column not found in ",lab," summary statistics",sep=""))
    }
    XtY=ss$Coef
    outnm="Coef"
  }else{
    stop('The data type indicator (Tind or Pind) need to be one of "GWAS", "Corr", or "Coef"')
  }
  
  # remove NA
  pos_na=is.na(XtY)
  XtY=XtY[!pos_na]
  ss=ss[!pos_na,]
  
  # subset markers based on reference
  nm=paste(ss$CHR,ss$BP,ss$A1,ss$A2,sep="_")
  nm_flp=paste(ss$CHR,ss$BP,ss$A2,ss$A1,sep="_")
  New_ss=data.frame(KG[,2], KG[,1], KG[,4], KG[,5], KG[,6],0)
  colnames(New_ss)=c("nm","CHR","BP","A1","A2",outnm)
  
  pos=New_ss$nm%in%nm
  pos_flp=New_ss$nm%in%nm_flp
  New_ss[pos,6]=XtY[match(New_ss[pos,1],nm)]
  New_ss[pos_flp,6]=-XtY[match(New_ss[pos_flp,1],nm_flp)]
  
  writeLines(strwrap(paste(sum(pos)+sum(pos_flp)," SNPs passed filtering for ",lab,sep="")))
  
  return(New_ss)
}
