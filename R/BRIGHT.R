Sol_path.plot <- function(out,highlight=NA,eta.plot=NA,lambda.plot=NA){
  
  non_zero_pos=highlight
  
  eta_vec=out[["eta_vec"]]
  lambda=out[["lambda"]]
  p=dim(out[["0"]]$Beta)[1]
  
  if(!is.na(eta.plot)){
    Beta=out[[as.character(eta.plot)]]$Beta
    cut_lambda_ind=which(colSums(Beta!=0)!=0)
    cut_Beta=Beta[,cut_lambda_ind]
    cut_lambda=lambda[cut_lambda_ind]
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(-log(max(cut_lambda)),-log(min(cut_lambda))),xlab = expression(-log(lambda)),ylab = "Estimated effect sizes",main=paste("lambda solution path (eta=",round(eta.plot,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }
  
  if(!is.na(lambda.plot)){
    cut_eta=eta_vec[-22]
    
    rw=which(lambda==lambda.plot)
    
    cut_Beta=c()
    for(i in eta_vec[-22]){
      cut_Beta=cbind(cut_Beta,out[[as.character(i)]][,rw])
    }
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(min(cut_eta),max(cut_eta)),xlab = expression(eta),ylab = "Estimated effect sizes", main=paste("eta solution path (lambda=",round(lambda.plot,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }
}

Sol_path.plot <- function(out,highlight=NA,eta.plot=NA,lambda.plot=NA){
  
  non_zero_pos=highlight
  
  eta_vec=out[["eta_vec"]]
  lambda=out[["lambda"]]
  p=dim(out[["0"]]$Beta)[1]
  
  if(!is.na(eta.plot)){
    Beta=out[[as.character(eta.plot)]]$Beta
    cut_lambda_ind=which(colSums(Beta!=0)!=0)
    cut_Beta=Beta[,cut_lambda_ind]
    cut_lambda=lambda[cut_lambda_ind]
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(-log(max(cut_lambda)),-log(min(cut_lambda))),xlab = expression(-log(lambda)),ylab = "Estimated effect sizes",main=paste("lambda solution path (eta=",round(eta.plot,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(-log(cut_lambda),cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }
  
  if(!is.na(lambda.plot)){
    cut_eta=eta_vec[-22]
    
    rw=which.min(abs(lambda-lambda.plot))[1]
    
    cut_Beta=c()
    for(i in eta_vec[-22]){
      cut_Beta=cbind(cut_Beta,out[[as.character(i)]]$Beta[,rw])
    }
    
    par(mar = c(5, 5, 3, 0.5))
    plot(1,type = "n",ylim = c(min(cut_Beta),max(cut_Beta)),xlim = c(min(cut_eta),max(cut_eta)),xlab = expression(eta),ylab = "Estimated effect sizes", main=paste("eta solution path (lambda=",round(lambda.plot,digits = 3),")",sep = ""), frame.plot = FALSE,cex.lab=1.5,cex.axis=1.5)
    for (i in 1:p) {
      if (!i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col=ggplot2::alpha("Grey",0.5),lty=1)}
    }
    for (i in 1:p) {
      if (i%in%non_zero_pos ){
        lines(cut_eta,cut_Beta[i,],ty="l",col="red",lwd=2)}
    }
  }
}