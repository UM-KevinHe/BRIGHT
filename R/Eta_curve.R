Eta_curve.plot <- function(Val,valid.type="AIC"){
  if(valid.type=="AIC"){
    eta_vec=Val[["eta_vec"]]
    lambda=log(Val[["lambda"]])
    AIC=Val[["AIC"]]
    dinom=dim(AIC)[1]
    
    m=matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
    
    layout(mat = m,heights = c(0.1,0.9))
    
    par(mar = c(0, 0, 0, 0))
    plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
    legend(x="center", legend=c("BRIGHT","Local"),col=c("black", "black"),lty=c(1,NA),lwd=c(3,NA),pch = c(16,15),cex = c(1.5,1.5),bty = "n",horiz = T)
    
    BRCor=AIC
    eta=eta_vec
    
    ind=which.min(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[rw,]),max(ER[rw,])),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = "APD",cex.lab=1.5,cex.axis=1.5)
    lines(eta,ER[rw,],lwd=3)
    points(eta[-1],ER[rw,-1],pch=16,cex=2)
    points(eta[1],ER[rw,1],pch=15,cex=2)
    
    cut_lambda=lambda
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[cut_lambda,cl]),max(ER[cut_lambda,cl])),xlim = c(min(-lambda[cut_lambda]),max(-lambda[cut_lambda])),xlab = expression(-log(lambda)),ylab = "APD",cex.lab=1.5,cex.axis=1.5)
    lines(-lambda[cut_lambda],ER[cut_lambda,cl],lwd=3)
    points(-lambda[cut_lambda],ER[cut_lambda,cl],pch=16,cex=2)
    
  }else if(valid.type=="Ind"){
    eta_vec=Val[["eta_vec"]]
    lambda=log(Val[["lambda"]])
    dinom=dim(Val[["BRMSE_rst"]])[1]
    
    m=matrix(c(1,1,2,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)
    
    layout(mat = m,heights = c(0.1,0.45,0.45))
    
    par(mar = c(0, 0, 0, 0))
    plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
    legend(x="center", legend=c("BRIGHT","Local","Prior"),col=c("black", "black","black"),lty=c(1,NA,NA),lwd=c(3,NA,NA),pch = c(16,15,17),cex = c(1.5,1.5,1.5),bty = "n",horiz = T)
    
    BRCor=Val[["BRMSE_rst"]]
    eta=eta_vec
    
    ind=which.min(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor
    
    ER_Exom=Val[["LASMSE"]]
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[rw,]),max(ER[rw,],ER_Exom)),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = "PD",cex.lab=1.5,cex.axis=1.5)
    lines(eta,ER[rw,],lwd=3)
    points(eta[-1],ER[rw,-1],pch=16,cex=2)
    points(eta[1],ER[rw,1],pch=15,cex=2)
    points(eta[21],ER_Exom,pch=17,cex=2)
    
    #cut_lambda=which(colSums(Val[["Beta_MSE"]]!=0)==0)[1]-1
    cut_lambda=length(lambda)
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[1:cut_lambda,cl]),max(ER[1:cut_lambda,cl],ER_Exom)),xlim = c(min(-lambda[1:cut_lambda]),max(-lambda[1:cut_lambda])),xlab = expression(-log(lambda)),ylab = "PD",cex.lab=1.5,cex.axis=1.5)
    lines(-lambda[1:cut_lambda],ER[1:cut_lambda,cl],lwd=3)
    points(-lambda[1:cut_lambda],ER[1:cut_lambda,cl],pch=16,cex=2)
    
    BRCor=Val[["BRCor_rst"]]
    eta=eta_vec
    
    ind=which.max(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor^2
    
    ER_Exom=Val[["LAScor"]]^2
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[rw,],ER_Exom),max(ER[rw,],ER_Exom)),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = expression(R^2),cex.lab=1.5,cex.axis=1.5)
    lines(eta,ER[rw,],lwd=3)
    points(eta[-1],ER[rw,-1],pch=16,cex=2)
    points(eta[1],ER[rw,1],pch=15,cex=2)
    points(eta[length(eta)],ER_Exom,pch=17,cex=2)
    
    #cut_lambda=which(colSums(Val[["Beta_Cor"]]!=0)==0)[1]-1
    cut_lambda=length(lambda)
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[1:cut_lambda,cl]),max(ER[1:cut_lambda,cl],ER_Exom)),xlim = c(min(-lambda[1:cut_lambda]),max(-lambda[1:cut_lambda])),xlab = expression(-log(lambda)),ylab = "R2",cex.lab=1.5,cex.axis=1.5)
    lines(-lambda[1:cut_lambda],ER[1:cut_lambda,cl],lwd=3)
    points(-lambda[1:cut_lambda],ER[1:cut_lambda,cl],pch=16,cex=2)
  }else if(valid.type=="Sum"){
    eta_vec=Val[["eta_vec"]]
    lambda=log(Val[["lambda"]])
    dinom=dim(Val[["BRMSE_rst"]])[1]
    
    m=matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
    
    layout(mat = m,heights = c(0.1,0.9))
    
    par(mar = c(0, 0, 0, 0))
    plot(1,type = "n",axes = F,xlab = "",ylab = "",frame.plot = FALSE)
    legend(x="center", legend=c("BRIGHT","Local","Prior"),col=c("black", "black","black"),lty=c(1,NA,NA),lwd=c(3,NA,NA),pch = c(16,15,17),cex = c(1.5,1.5,1.5),bty = "n",horiz = T)
    
    BRCor=Val[["BRMSE_rst"]]
    eta=eta_vec
    
    ind=which.min(BRCor)
    
    rw=ind%%dinom
    cl=ind%/%dinom+1
    ER=BRCor
    
    ER_Exom=Val[["LASMSE"]]
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[rw,]),max(ER[rw,],ER_Exom)),xlim = c(min(eta),max(eta)),xlab = expression(eta),ylab = "APD",cex.lab=1.5,cex.axis=1.5)
    lines(eta,ER[rw,],lwd=3)
    points(eta[-1],ER[rw,-1],pch=16,cex=2)
    points(eta[1],ER[rw,1],pch=15,cex=2)
    points(eta[length(eta)],ER_Exom,pch=17,cex=2)
    
    cut_lambda=which(colSums(Val[["Beta_MSE"]]!=0)!=0)
    
    par(mar = c(5, 5, 0.5, 0.5))
    plot(1,type = "n",ylim = c(min(ER[cut_lambda,cl]),max(ER[cut_lambda,cl],ER_Exom)),xlim = c(min(-lambda[cut_lambda]),max(-lambda[cut_lambda])),xlab = expression(-log(lambda)),ylab = "APD",cex.lab=1.5,cex.axis=1.5)
    lines(-lambda[cut_lambda],ER[cut_lambda,cl],lwd=3)
    points(-lambda[cut_lambda],ER[cut_lambda,cl],pch=16,cex=2)
  }else{
    stop('valid.type should be either "AIC", "Ind", or "Sum" representing the AIC criterion, individual-level data, or summary-level data')
  }
}