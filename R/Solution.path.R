Solution.path <- function(Val,out,highlight=NA,criteria="R2"){
  non_zero_pos=highlight
  
  eta_vec=Val[["eta_vec"]]
  lambda=Val[["lambda"]]
  dinom=dim(Val[["BRMSE_rst"]])[1]
  p=dim(out[[as.character(Val[["Best_eta_MSE"]])]]$Beta)[1]
  
  m=matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
  
  layout(mat = m,heights = c(0.1,0.9))
  
  par(mar = c(0, 0, 0, 0))
  plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
  legend(x="center", legend=c("Highlight markers","Other markers"),col=c("red", ggplot2::alpha("Grey",0.5)),lty=c(1,1),lwd=c(2,1),cex = c(1.5,1.5),bty = "n",horiz = T)
  
  if(criteria=="PD" || criteria=="APD"){
    Best_lambda=Val[["Best_lambda_MSE"]]
    Best_eta=Val[["Best_eta_MSE"]]
    
    Beta=out[[as.character(Val[["Best_eta_MSE"]])]]$Beta
    cut_lambda_ind=which(colSums(Beta!=0)!=0)
    cut_Beta=Beta[,cut_lambda_ind]
    cut_lambda=lambda[cut_lambda_ind]
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(-log(max(cut_lambda)),-log(min(cut_lambda))),xlab = expression(-log(lambda)),ylab = "Estimated effect sizes", main=paste("lambda solution path (eta=",round(Best_eta,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
    
    BRCor=Val[["BRMSE_rst"]]
    cut_eta=eta_vec
    
    ind=which.min(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor
    
    ER_Exom=Val[["LASMSE"]]
    
    cut_Beta=c()
    for(i in eta_vec){
      cut_Beta=cbind(cut_Beta,out[[as.character(i)]]$Beta[,rw])
    }
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(min(cut_eta),max(cut_eta)),xlab = expression(eta),ylab = "Estimated effect sizes", main=paste("eta solution path (lambda=",round(Best_lambda,digits = 3),")",sep = ""),frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }else if(criteria=="R2"){
    print("***")
    Best_lambda=Val[["Best_lambda_Cor"]]
    Best_eta=Val[["Best_eta_Cor"]]
    
    Beta=out[[as.character(Val[["Best_eta_Cor"]])]]$Beta
    cut_lambda_ind=which(colSums(Beta!=0)!=0)
    cut_Beta=Beta[,cut_lambda_ind]
    cut_lambda=lambda[cut_lambda_ind]
    
    non_zero_pos=highlight
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(-log(max(cut_lambda)),-log(min(cut_lambda))),xlab = expression(-log(lambda)),ylab = "Estimated effect sizes", main=paste("lambda solution path (eta=",round(Best_eta,digits = 3),")",sep = ""),frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
    
    BRCor=Val[["BRCor_rst"]]
    cut_eta=eta_vec
    
    ind=which.max(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor^2
    
    ER_Exom=Val[["LAScor"]]^2
    
    cut_Beta=c()
    for(i in eta_vec){
      cut_Beta=cbind(cut_Beta,out[[as.character(i)]]$Beta[,rw])
    }
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(min(cut_eta),max(cut_eta)),xlab = expression(eta),ylab = "Estimated effect sizes", main=paste("eta solution path (lambda=",round(Best_lambda,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }else{
    stop('criteria must be either "PD", "APD", or "R2" representing predictive deviance, approximated predictive deviance, or R2')
  }
  
}
