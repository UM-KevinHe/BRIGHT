Valid.Ind <- function(out, Testpheno, Testgeno){
  gene_bed = BEDMatrix::BEDMatrix(paste(Testgeno,".bed",sep=""))
  gene_bed=as.matrix(gene_bed)
  XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
  XB_Exom=XG$X%*%out[["Prior"]]
  phe=Testpheno[,3]
  phe_std=(phe-mean(phe))/sd(phe)
  LAScor=cor(XB_Exom, phe_std)
  LASMSE=colMeans((XB_Exom-as.vector(phe_std))^2)
  eta_vec=out[["eta_vec"]]
  AIC=out$AIC
  
  value=rep(0,6)
  pos=rep(0,6)
  eta_best=rep(0,3)
  ER=c()
  BRMSE_rst=c()
  BRCor_rst=c()
  Vali_MSE_rst=c()
  Vali_cor_rst=c()
  Beta=list()
  Vali_MSE_best=Inf
  Vali_cor_best=-Inf
  AIC_best=which.min(AIC)
  XB_MSE=c()
  XB_Cor=c()
  Beta_MSE=c()
  Beta_Cor=c()
  
  phe_test_std=(phe-mean(phe))/sd(phe)
  beta=out[[as.character(0)]]$Beta
  XB_Ind2=XG$X%*%beta
  Vali_MSE=colMeans((XB_Ind2-as.vector(phe_test_std))^2)
  Vali_cor=cor(XB_Ind2, phe_test_std)
  XB_local_Cor=XB_Ind2[,which.max(Vali_cor)]
  Beta_local_Cor=beta[,which.max(Vali_cor)]
  XB_local_MSE=XB_Ind2[,which.min(Vali_MSE)]
  Beta_local_MSE=beta[,which.min(Vali_MSE)]
  
  for(eta in eta_vec){
    beta=out[[as.character(eta)]]$Beta
    
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
      Beta_MSE=beta
      Beta_BRIGHT_MSE=beta[,which.min(Vali_MSE)]
    }
    if(max(Vali_cor,na.rm=T)>Vali_cor_best){
      Vali_cor_best=max(Vali_cor,na.rm=T)
      value[3]=BRcor[Vali_cor_ind]
      value[4]=BRMSE[Vali_cor_ind]
      eta_best[2]=eta
      XB_Cor=XB_Ind2[,which.max(Vali_cor)]
      Beta_Cor=beta
      Beta_BRIGHT_Cor=beta[,which.max(Vali_cor)]
    }
    
    Vali_MSE_rst=cbind(Vali_MSE_rst,Vali_MSE)
    Vali_cor_rst=cbind(Vali_cor_rst,Vali_cor)
    
    BRMSE_rst=cbind(BRMSE_rst,BRMSE)
    BRCor_rst=cbind(BRCor_rst,BRcor)
    
  }
  
  dinom=dim(BRCor_rst)[1]
  ind=which.max(BRCor_rst)
  rw=ind%%dinom
  cl=ind%/%dinom+1
  if(rw==0){
    rw=dinom
    cl=cl-1
  }
  Best_eta_Cor=eta_vec[cl]
  Best_lambda_Cor=out[["lambda"]][rw]
  
  dinom=dim(BRMSE_rst)[1]
  ind=which.min(BRMSE_rst)
  rw=ind%%dinom
  cl=ind%/%dinom+1
  if(rw==0){
    rw=dinom
    cl=cl-1
  }
  Best_eta_MSE=eta_vec[cl]
  Best_lambda_MSE=out[["lambda"]][rw]
  
  Best_eta_AIC=out[["Best_eta_AIC"]]
  Best_lambda_AIC=out[["Best_lambda_AIC"]]
  Beta_BRIGHT_AIC=out[["Beta_BRIGHT_AIC"]]
  Beta_local_AIC=out[["Beta_local_AIC"]]
  
  print(paste("Best eta based on R2 is",round(Best_eta_Cor,digits = 3)))
  print(paste("Best lambda based on R2 is",round(Best_lambda_Cor,digits = 3)))
  print(paste("Best eta based on MSPE is",round(Best_eta_MSE,digits = 3)))
  print(paste("Best lambda based on MSPE is",round(Best_lambda_MSE,digits = 3)))
  print(paste("Best eta based on AIC is",round(Best_eta_AIC,digits = 3)))
  print(paste("Best lambda based on AIC is",round(Best_lambda_AIC,digits = 3)))
  
  rst=list()
  rst[["BRMSE_rst"]]=BRMSE_rst
  rst[["BRCor_rst"]]=BRCor_rst
  rst[["AIC"]]=AIC
  rst[["LAScor"]]=LAScor
  rst[["LASMSE"]]=LASMSE
  rst[["XB_Exom"]]=XB_Exom
  rst[["Beta_Exom"]]=out[["Prior"]]
  rst[["XB_local_MSE"]]=XB_local_MSE
  rst[["XB_local_Cor"]]=XB_local_Cor
  rst[["Beta_local_MSE"]]=Beta_local_MSE
  rst[["Beta_local_Cor"]]=Beta_local_Cor
  rst[["XB_MSE"]]=XB_MSE
  rst[["XB_Cor"]]=XB_Cor
  rst[["Beta_MSE"]]=Beta_MSE
  rst[["Beta_Cor"]]=Beta_Cor
  rst[["Beta_BRIGHT_MSE"]]=Beta_BRIGHT_MSE
  rst[["Beta_BRIGHT_Cor"]]=Beta_BRIGHT_Cor
  rst[["Beta_BRIGHT_AIC"]]=Beta_BRIGHT_AIC
  rst[["Beta_local_AIC"]]=Beta_local_AIC
  rst[["phe"]]=phe
  rst[["phe_std"]]=phe_std
  rst[["eta_vec"]]=eta_vec
  rst[["lambda"]]=out[["lambda"]]
  rst[["Best_eta_Cor"]]=Best_eta_Cor
  rst[["Best_lambda_Cor"]]=Best_lambda_Cor
  rst[["Best_eta_MSE"]]=Best_eta_MSE
  rst[["Best_lambda_MSE"]]=Best_lambda_MSE
  rst[["Best_eta_AIC"]]=Best_eta_AIC
  rst[["Best_lambda_AIC"]]=Best_lambda_AIC
  return(rst)
}
