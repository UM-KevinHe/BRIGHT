Valid.Sum <- function(out, ValiCorr, ValiRef){
  gene_bed = BEDMatrix::BEDMatrix(paste(ValiRef,".bed",sep=""))
  gene_bed=as.matrix(gene_bed)
  XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
  XB_Exom=XG$X%*%out[["Prior"]]
  LASMSE=t(XB_Exom)%*%XB_Exom/nrow(XG$X)-2*ValiCorr$Corr%*%out[["Prior"]]
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
  Beta_MSE=c()
  Beta_Cor=c()
  
  beta=out[[as.character(0)]]$Beta
  XB_Ind2=XG$X%*%beta
  Vali_MSE=diag(t(XB_Ind2)%*%XB_Ind2/nrow(XG$X))-2*ValiCorr$Corr%*%beta
  XB_local_MSE=XB_Ind2[,which.min(Vali_MSE)]
  Beta_local_MSE=beta[,which.min(Vali_MSE)]
  
  for(eta in eta_vec){
    beta=out[[as.character(eta)]]$Beta
    
    XB_Ind2=XG$X%*%beta
    
    Vali_MSE=diag(t(XB_Ind2)%*%XB_Ind2/nrow(XG$X))-2*ValiCorr$Corr%*%beta
    Vali_MSE_ind=which.min(Vali_MSE)
    
    BRMSE=diag(t(XB_Ind2)%*%XB_Ind2/nrow(XG$X))-2*ValiCorr$Corr%*%beta
    
    if(min(Vali_MSE,na.rm=T)<Vali_MSE_best){
      Vali_MSE_best=min(Vali_MSE,na.rm=T)
      value[2]=BRMSE[Vali_MSE_ind]
      eta_best[1]=eta
      XB_MSE=XB_Ind2[,which.min(Vali_MSE)]
      Beta_MSE=beta
    }
    
    Vali_MSE_rst=cbind(Vali_MSE_rst,Vali_MSE)
    
    BRMSE_rst=cbind(BRMSE_rst,t(BRMSE))
    
  }
  
  dinom=dim(BRMSE_rst)[1]
  ind=which.min(BRMSE_rst)
  rw=ind%%dinom
  cl=ind%/%dinom+1
  Best_eta_MSE=eta_vec[cl]
  Best_lambda_MSE=out[["lambda"]][rw]
  
  Best_eta_AIC=out[["Best_eta_AIC"]]
  Best_lambda_AIC=out[["Best_lambda_AIC"]]
  Beta_BRIGHT_AIC=out[["Beta_BRIGHT_AIC"]]
  Beta_local_AIC=out[["Beta_local_AIC"]]
  
  print(paste("Best eta based on AMSPE is",round(Best_eta_MSE,digits = 3)))
  print(paste("Best lambda based on AMSPE is",round(Best_lambda_MSE,digits = 3)))
  print(paste("Best eta based on AIC is",round(Best_eta_AIC,digits = 3)))
  print(paste("Best lambda based on AIC is",round(Best_lambda_AIC,digits = 3)))
  
  
  rst=list()
  rst[["BRMSE_rst"]]=BRMSE_rst
  rst[["LASMSE"]]=LASMSE
  rst[["Beta_Exom"]]=out[["Prior"]]
  rst[["Beta_local_MSE"]]=Beta_local_MSE
  rst[["Beta_MSE"]]=Beta_MSE
  rst[["Beta_BRIGHT_AIC"]]=Beta_BRIGHT_AIC
  rst[["Beta_local_AIC"]]=Beta_local_AIC
  rst[["eta_vec"]]=eta_vec
  rst[["lambda"]]=out[["lambda"]]
  rst[["Best_eta_MSE"]]=Best_eta_MSE
  rst[["Best_lambda_MSE"]]=Best_lambda_MSE
  rst[["Best_eta_AIC"]]=Best_eta_AIC
  rst[["Best_lambda_AIC"]]=Best_lambda_AIC
  return(rst)
}