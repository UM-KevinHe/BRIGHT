ROC.plot <- function(Val, Pct=0.5, criteria="R2"){
  if(criteria=="R2"){
    if(is.null(Val[["XB_local_Cor"]])){
      warning('"R2" not found using "AIC" instead')
      XB_local=Val[["XB_local_AIC"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_AIC"]]
    }else{
      XB_local=Val[["XB_local_Cor"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_Cor"]]
    }
    
  }else if(criteria=="PD"){
    if(is.null(Val[["XB_local_Cor"]])){
      warning('"PD" not found using "AIC" instead')
      XB_local=Val[["XB_local_AIC"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_AIC"]]
    }else{
      XB_local=Val[["XB_local_MSE"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_MSE"]]
    }
    
  }else if(criteria=="APD"){
    if(is.null(Val[["XB_local_Cor"]])){
      warning('"APD" not found using "AIC" instead')
      XB_local=Val[["XB_local_AIC"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_AIC"]]
    }else{
      XB_local=Val[["XB_local_MSE"]]
      XB_prior=Val[["XB_Exom"]]
      XB_BRIGHT=Val[["XB_MSE"]]
    }
    
  }else if(criteria=="AIC"){
    XB_local=Val[["XB_local_AIC"]]
    XB_prior=Val[["XB_Exom"]]
    XB_BRIGHT=Val[["XB_AIC"]]
  }else{
    stop('criteria must be either "AIC", "PD", "APD", or "R2" representing AIC criterion, predictive deviance, approximated predictive deviance, or R2')
  }
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
  legend("bottomright", legend=c(paste("BRIGHT (AUC:",round(as.numeric(ROC_BRIGHT$auc),digits = 3),")",sep = ""),paste("Local (AUC:",round(as.numeric(ROC_local$auc),digits = 3),")",sep = ""),paste("Prior (AUC:",round(as.numeric(ROC_prior$auc),digits = 3),")",sep = "")),col=c("black", "blue", "red"),lty=c(1,3,2), lwd = c(3,2,2), cex=1)
}