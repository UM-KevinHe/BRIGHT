Density.plot <- function(Val,Pct,criteria="R2"){
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
  
  ind_local=XB_local>quantile(XB_local,probs = c(0,0.1,0.5,Pct,1))[4]
  ind_prior=XB_prior>quantile(XB_prior,probs = c(0,0.1,0.5,Pct,1))[4]
  ind_BRIGHT=XB_BRIGHT>quantile(XB_BRIGHT,probs = c(0,0.1,0.5,Pct,1))[4]
  
  high_local=phe[ind_local]
  low_local=phe[!ind_local]
  
  high_prior=phe[ind_prior]
  low_prior=phe[!ind_prior]
  
  high_BRIGHT=phe[ind_BRIGHT]
  low_BRIGHT=phe[!ind_BRIGHT]
  
  m=matrix(c(1,1,1,2,3,4),nrow = 2,ncol = 3,byrow = TRUE)
  
  layout(mat = m,heights = c(0.1,0.9))
  
  par(mar = c(0, 0, 0, 0))
  plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
  legend(x="center", legend=c(paste("Upper ",round(Pct*100)," PRS percentile",sep = ""),paste("Lower ",round(Pct*100)," PRS percentile",sep = "")),col=c("red", "black"),lty=c(2,1),lwd=c(3,3),cex = c(2,2),bty = "n",horiz = T)
  
  dist_local=density(high_local)$x[which.max(density(high_local)$y)]-density(low_local)$x[which.max(density(low_local)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.2),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("Local (",round(dist_local,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_local),lty=2,lwd=3,col="red")
  lines(density(low_local),lty=1,lwd=3)
  abline(v = density(high_local)$x[which.max(density(high_local)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_local)$x[which.max(density(low_local)$y)], col="black", lwd=3, lty=3)
  
  dist_prior=density(high_prior)$x[which.max(density(high_prior)$y)]-density(low_prior)$x[which.max(density(low_prior)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.2),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("Prior (",round(dist_prior,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_prior),lty=2,lwd=3,col="red")
  lines(density(low_prior),lty=1,lwd=3)
  abline(v = density(high_prior)$x[which.max(density(high_prior)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_prior)$x[which.max(density(low_prior)$y)], col="black", lwd=3, lty=3)
  
  dist_BRIGHT=density(high_BRIGHT)$x[which.max(density(high_BRIGHT)$y)]-density(low_BRIGHT)$x[which.max(density(low_BRIGHT)$y)]
  par(mar = c(5, 5, 3, 0.5))
  plot(1,type = "n",ylim = c(0,0.2),xlim = c(min(phe),max(phe)),xlab = "Outcome",ylab = "Density",cex.lab=2,cex.axis=2,main=paste("BRIGHT (",round(dist_BRIGHT,digits = 2),")",sep = ""),cex.main=2,frame.plot=F)
  box(bty="l")
  lines(density(high_BRIGHT),lty=2,lwd=3,col="red")
  lines(density(low_BRIGHT),lty=1,lwd=3)
  abline(v = density(high_BRIGHT)$x[which.max(density(high_BRIGHT)$y)], col="red", lwd=3, lty=3)
  abline(v = density(low_BRIGHT)$x[which.max(density(low_BRIGHT)$y)], col="black", lwd=3, lty=3)
}