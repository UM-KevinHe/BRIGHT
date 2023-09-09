Test.Ind <- function(Val, Testpheno, Testgeno){
  if(!is.null(Val[["Beta_local_Cor"]])){
    gene_bed = BEDMatrix::BEDMatrix(paste(Testgeno,".bed",sep=""))
    gene_bed=as.matrix(gene_bed)
    XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
    XB_Exom=XG$X%*%Val[["Beta_Exom"]]
    XB_local_MSE=XG$X%*%Val[["Beta_local_MSE"]]
    XB_local_Cor=XG$X%*%Val[["Beta_local_Cor"]]
    XB_MSE=XG$X%*%Val[["Beta_BRIGHT_MSE"]]
    XB_Cor=XG$X%*%Val[["Beta_BRIGHT_Cor"]]
    XB_local_AIC=XG$X%*%Val[["Beta_local_AIC"]]
    XB_AIC=XG$X%*%Val[["Beta_BRIGHT_AIC"]]
    phe=Testpheno[,3]
    
    rst=list()
    rst[["XB_Exom"]]=XB_Exom
    rst[["XB_local_MSE"]]=XB_local_MSE
    rst[["XB_local_Cor"]]=XB_local_Cor
    rst[["XB_MSE"]]=XB_MSE
    rst[["XB_Cor"]]=XB_Cor
    rst[["XB_local_AIC"]]=XB_local_AIC
    rst[["XB_AIC"]]=XB_AIC
    rst[["phe"]]=phe
    return(rst)
  }else if(!is.null(Val[["Beta_local_Cor"]])){
    gene_bed = BEDMatrix::BEDMatrix(paste(Testgeno,".bed",sep=""))
    gene_bed=as.matrix(gene_bed)
    XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
    XB_Exom=XG$X%*%Val[["Beta_Exom"]]
    XB_local_MSE=XG$X%*%Val[["Beta_local_MSE"]]
    XB_MSE=XG$X%*%Val[["Beta_BRIGHT_MSE"]]
    XB_local_AIC=XG$X%*%Val[["Beta_local_AIC"]]
    XB_AIC=XG$X%*%Val[["Beta_BRIGHT_AIC"]]
    phe=Testpheno[,3]
    
    rst=list()
    rst[["XB_Exom"]]=XB_Exom
    rst[["XB_local_MSE"]]=XB_local_MSE
    rst[["XB_MSE"]]=XB_MSE
    rst[["XB_local_AIC"]]=XB_local_AIC
    rst[["XB_AIC"]]=XB_AIC
    rst[["phe"]]=phe
    return(rst)
  }else{
    gene_bed = BEDMatrix::BEDMatrix(paste(Testgeno,".bed",sep=""))
    gene_bed=as.matrix(gene_bed)
    XG <- newXG(gene_bed, 1:ncol(gene_bed), rep(1,ncol(gene_bed)), 1, FALSE)
    XB_Exom=XG$X%*%Val[["Prior"]]
    XB_local_AIC=XG$X%*%Val[["Beta_local_AIC"]]
    XB_AIC=XG$X%*%Val[["Beta_BRIGHT_AIC"]]
    phe=Testpheno[,3]
    
    rst=list()
    rst[["XB_Exom"]]=XB_Exom
    rst[["XB_local_AIC"]]=XB_local_AIC
    rst[["XB_AIC"]]=XB_AIC
    rst[["phe"]]=phe
    return(rst)
  }
  
}
