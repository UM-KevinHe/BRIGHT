#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <string.h>
using namespace Rcpp;

//' Second Norm
//'
//' Calculate second norm (Euclidean distance) of a vertical vector.
//'
//'@param x a numeric column vector to calculate the second norm.
//'@param p an integer represents the dimension of x.
//'@return A scalar of the computed second norm of x.
//'@examples
//' x <- t(t(rnorm(100)))
//' Norm(x,length(x))
// [[Rcpp::export]]
double Norm(Eigen::MatrixXd x, int p) {
  double x_norm = 0;
  for (int j=0; j<p; j++) x_norm = x_norm + x(j,0)*x(j,0);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

//' Soft-thresholding
//'
//' Calculate the soft-thresholding result corresponding to the LASSO penalty.
//'
//'@param z a numeric scalar to be soft-thresholded.
//'@param l a numeric scalar for the soft-thresholding parameter.
//'@return If z <= l, 0 is outputted.
//'@return If z > l, z-l is outputted.
//'@examples 
//' S(5,1)
//' S(0.5,1)
// [[Rcpp::export]]
double S(double z, double l) {
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

//' Firm-thresholding
//'
//' Calculate the Firm-thresholding result corresponding to the MCP penalty.
//'
//'@param z a numeric scalar to be firm-thresholded.
//'@param l1 a numeric scalar for the first firm-thresholding parameter.
//'@param l2 a numeric scalar for the second firm-thresholding parameter.
//'@param gamma a numeric scalar for the third firm-thresholding parameter.
//'@return A firm-thresholded scalar of z.
//'@examples
//' Ff(5,1,1,1)
//' Ff(0.5,1,1,1)
// [[Rcpp::export]]
double Ff(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(1+l2-1/gamma));
  else return(z/(1+l2));
}

//' SCAD-thresholding
//'
//' Calculate the SCAD-thresholding result corresponding to the SCAD penalty.
//'
//'@param z a numeric scalar to be SCAD-thresholded.
//'@param l1 a numeric scalar for the first SCAD-thresholding parameter.
//'@param l2 a numeric scalar for the second SCAD-thresholding parameter.
//'@param gamma a numeric scalar for the third SCAD-thresholding parameter.
//'@return A SCAD-thresholded scalar of z.
//'@examples
//' Fs(5,1,1,1)
//' Fs(0.5,1,1,1)
// [[Rcpp::export]]
double Fs(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(1+l2));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(1-1/(gamma-1)+l2));
  else return(z/(1+l2));
}

//' LD Block Structured SNPs grouping
//'
//' Calculate the LD grouping of input SNPs based on existing block structures (Berisa et. al.). The chr and pos must be matched and sorted.
//'
//'@param chr a numeric verticle vector of the chromosomes each SNP belong to.
//'@param pos a numeric verticle vector of the position in base pair (BP) of each SNP.
//'@param LDB a numeric matrix for the block structure information directly read in from Berisa's bed file.
//'@return A numeric verticle vector, representing the block structured grouping of the current SNP list; SNP with (LD[i],LD[i+1]) are to be within one group.
//'@examples
//' chr=t(t(c(1,1,1,2,2,3,3)))
//' pos=t(t(c(10,20,30,1,20,9,20)))
//' LDB=data.frame(chr=c(1,1,1,2,2,2,3,3,3),start=c(1,15,40,1,30,60,1,50,100),stop=c(15,40,80,30,60,100,50,100,150))
//' LDB=as.matrix(LDB)
//' LD(chr,pos,LDB)
// [[Rcpp::export]]
Eigen::MatrixXd LD(Eigen::MatrixXd chr,Eigen::MatrixXd pos,Eigen::MatrixXd LDB){
  int j=-1;
  int p=pos.rows();
  int chr_mark=1;
  Eigen::MatrixXd blk(1,1);
  blk.setZero(1,1);
  while(1==1){
    j++;
    if(pos(0,0)<LDB(j,2)){
      break;
    }
  }

  for(int i=0;i<p;i++){
    if(chr(i,0)!=chr_mark){
      Eigen::MatrixXd blk_tmp=blk;
      blk.resize(blk.rows()+1,1);
      blk.block(0,0,blk.rows()-1,1)=blk_tmp;
      blk(blk.rows()-1,0)=i;
      chr_mark=chr(i,0);
      while(1==1){
        j++;
        if( (chr(i,0)==LDB(j,0)) and (pos(i,0)<LDB(j,2)) ){
          break;
        }
      }
    }
    if(pos(i,0)>=LDB(j,2)){
      Eigen::MatrixXd blk_tmp=blk;
      blk.resize(blk.rows()+1,1);
      blk.block(0,0,blk.rows()-1,1)=blk_tmp;
      blk(blk.rows()-1,0)=i;
      while(1==1){
        j++;
        if(pos(i,0)<LDB(j,2)){
          break;
        }
      }
    }
  }
  return(blk);
}

/*
// [[Rcpp::export]]
int Sigma_LD(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXd& blk, double tau, int P, int N){
  int lblk=blk.rows();
  for(int k=0;k<lblk-1;k++){
    for(int i=blk(k,0);i<blk(k+1,0);i++){
      Rcout << i;
      Rcout << "\n";
      Sig.insert(i,i) = 1.0;
      for(int j=i+1;j<blk(k+1,0);j++){
        double tmp=S((X.block(0,i,N,1).transpose()*X.block(0,j,N,1))(0,0)/N,tau);
        if(tmp!=0){
          Sig.insert(i,j) = tmp;
          Sig.insert(j,i) = tmp;
        }
      }
    }
  }
  Sig.makeCompressed();
  return(0);
}

// [[Rcpp::export]]
int Sigma_LD(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXd& blk, double tau, int P, int N){
  int lblk=blk.rows();
  for(int k=0;k<lblk-1;k++){
    Eigen::MatrixXd mat=(X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)).transpose()*X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)))/N;
    for(int i=blk(k,0);i<blk(k+1,0);i++){
      Rcout << i;
      Rcout << "\n";
      Sig.insert(i,i) = 1.0;
      for(int j=i+1;j<blk(k+1,0);j++){
        double tmp=S(mat(i-blk(k,0),j-blk(k,0)),tau);
        if(tmp!=0){
          Sig.insert(i,j) = tmp;
          Sig.insert(j,i) = tmp;
        }
      }
    }
  }
  Sig.makeCompressed();
  return(0);
}

// [[Rcpp::export]]
int Sigma_LD(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXi& blk, double tau, int P, int N){
  int lblk=blk.rows();
  for(int k=0;k<lblk-1;k++){
    Eigen::MatrixXd mat(P,blk(k+1,0)-blk(k,0));
    //mat.setZero(P,blk(k+1,0)-blk(k,0));
    //mat.block(0,0,blk(k+1,0)-blk(k,0),blk(k+1,0)-blk(k,0))=(X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)).transpose()*X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)))/N;
    Sig.middleCols(blk(k,0),blk(k+1,0)-blk(k,0))=mat.block(0,0,P,blk(k+1,0)-blk(k,0));
  }
  Sig.makeCompressed();
  return(0);
}
*/

//' Block Structured LD Calculation
//'
//' Calculate the LD based on existing block structures.
//'
//'@param X a numeric matrix of the standardized genotypes outputted through newXG() function.
//'@param Sig a numeric sparse matrix pointer, indicating the memory location where the output sparse LD matrix should be stored.
//'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
//'@param tau a numeric scalar represents the element-wise soft-thresholding parameter for generating the LD.
//'@param P an integer scalar for the column dimension (genotype dimension) of X.
//'@param N an integer scalar for the row dimension (sample size) of X.
//'@return Write the calculated sparse LD matrix to the memory location provided in Sig.
// [[Rcpp::export]]
Eigen::SparseMatrix<double> Sigma_LD(Eigen::MatrixXd& X, Eigen::MatrixXi& blk, double tau, int P, int N){
  int lblk=blk.rows();
  
  Eigen::SparseMatrix<double, 0> Sig(P,P);
  Sig.reserve(Eigen::VectorXi::Constant(P,3000));
  
  typedef Eigen::Triplet<double> T;
  std::vector< T > TL;
  TL.reserve(10000000);
  
  for(int k=0;k<lblk-1;k++){
    Eigen::MatrixXd mat=(X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)).transpose()*X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)))/N;
      for(int i=blk(k,0);i<blk(k+1,0);i++){
      //Rcout << i;
      //Rcout << "\n";
      TL.push_back(T(i,i,1.0));
      for(int j=i+1;j<blk(k+1,0);j++){
        double tmp=S(mat(i-blk(k,0),j-blk(k,0)),tau);
        if(tmp!=0){
          TL.push_back(T(i,j,tmp));
          TL.push_back(T(j,i,tmp));
        }
      }
    }
  }
  Sig.setFromTriplets(TL.begin(), TL.end());
  return(Sig);
}


/*
// [[Rcpp::export]]
int Sigma_ST(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXi& blk, double tau, int P, int N){
  for(int i=0;i<P;i++){
    Rcout << i;
    Rcout << "\n";
    Sig.insert(i,i) = 1.0;
    for(int j=i+1;j<std::min(P,P);j++){
      double tmp=S((X.block(0,i,N,1).transpose()*X.block(0,j,N,1))(0,0)/N,tau);
      if(tmp!=0){
        Sig.insert(i,j) = tmp;
        Sig.insert(j,i) = tmp;
      }
    }
  }
  Sig.makeCompressed(); 
  return(0);
}
*/
/*
// [[Rcpp::export]]
int Sigma_ST(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXi& blk, double tau, int P, int N){
  typedef Eigen::Triplet<double> T;
  std::vector< T > TL;
  TL.reserve(10000000);

  for(int i=0;i<P;i++){
    Rcout << i;
    Rcout << "\n";
    TL.push_back(T(i,i,1.0));
    Eigen::MatrixXd mat(1,P);
    mat.setZero(1,P);
    mat.block(0,i+1,1,P-i-1)=(X.block(0,i,N,1).transpose()*X.block(0,i+1,N,P-i-1))/N;
    for(int j=i+1;j<P;j++){
      double tmp=S(mat(0,j),tau);
      if(tmp!=0){
          TL.push_back(T(i,j,tmp));
          TL.push_back(T(j,i,tmp));
      }
    }
  }
  Sig.setFromTriplets(TL.begin(), TL.end()); 
  return(0);
}
*/

// [[Rcpp::export]]
void gd_gaussian(Eigen::MatrixXd & a, Eigen::MatrixXd & b, Eigen::SparseMatrix<double> & Sig, Eigen::MatrixXd & K1,
                 Eigen::MatrixXd & Xtbt, Eigen::MatrixXd & XtY, Eigen::MatrixXd & df,
                 int & l, int & P, int & g, int & penalty, double & lam1, double & lam2,
                 double & gamma, double & eta, double & maxChange) {
  
  int Kg = K1(g+1,0) - K1(g,0);
  // Calculate z
  Eigen::MatrixXd z(Kg,1);
  z.setZero(Kg,1);
  //for (int j=K1(g,0); j<K1(g+1,0); j++) {
  //  z(j-K1(g,0),0) = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b.block(0,l,P,1))(0,0)+b(j,l);

  /*
  double tmp=(Sig.col(K1(g,0)).transpose()*b.block(0,l,P,1))(0,0)+b(K1(g,0),l);
  if (NumericVector::is_na(tmp)){
      Rcout << j;
      Rcout << "j\n";
      Rcout << "first part\n";
      Rcout << Sig.col(j).transpose();
      Rcout << "\nfirst part end\n";
      Rcout << "Second part\n";
      Rcout << b.block(0,l,P,1).transpose();
      Rcout << "\nSecond part end\n";
      Rcout << "Third part\n";
      Rcout << b(j,l);
      Rcout << "\nThird part end\n";
  }
  */
    
  //}

  z = (XtY.block(K1(g,0),0,Kg,1)+eta*Xtbt.block(K1(g,0),0,Kg,1))/(1+eta)-(Sig.middleCols(K1(g,0),Kg).transpose()*b.block(0,l,P,1))+b.block(K1(g,0),l,Kg,1);

  double z_norm = Norm(z, Kg);
  
  if (NumericVector::is_na(z_norm)){
      Rcout << K1(g,0);
      Rcout << "K1(g,0)\n";
  }


  // Update b
  double len;
  if (penalty==1) len = S(z_norm, lam1) / (1+lam2);
  if (penalty==2) len = Ff(z_norm, lam1, lam2, gamma);
  if (penalty==3) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 || a(K1(g,0),0) != 0) {
    // If necessary, update beta and r
    for (int j=K1(g,0); j<K1(g+1,0); j++) {
      b(j,l) = len * z(j-K1(g,0),0) / z_norm;
      double shift = b(j,l)-a(j,0);
      /*Rcout << "shift\n";
      Rcout << shift;
      Rcout << "shift_end\n";
      if(fabs(shift)>10){
      Rcout << "len\n";
      Rcout << len;
      Rcout << "len_end\n";
      Rcout << "z\n";
      Rcout << z(j-K1(g,0),0);
      Rcout << "z_end\n";
      Rcout << "z_norm\n";
      Rcout << z_norm;
      Rcout << "z_norm_end\n";
}*/

      if (fabs(shift) > maxChange) maxChange = fabs(shift);
    }
  }
  // Update df
  if (len > 0) df(l,0) += Kg * len / z_norm;

}

/*
// [[Rcpp::export]]
Eigen::SparseMatrix<double> Sigma_test(Eigen::MatrixXd& X, double tau, int P, int N){
  Eigen::SparseMatrix<double> Sig(P,P);
  Sig.reserve(Eigen::VectorXi::Constant(P,1000));
  Sigma(X, Sig, tau, P, N);
  return(Sig.col(1));
}
*/

//' Maximum Lambda selection
//'
//' Calculate the maximum lambda, which corresponding to penalizing all coefficients to be zero
//'
//'@param XtY a numeric verticle vector of the marginal correlation between genotype and outcome, can be recovered from GWAS summary statistics by function p2cor().
//'@param tilde_beta a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
//'@param X a numeric matrix of reference genotypes, default is from 1000 genome project.
//'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
//'@param K1 a numeric verticle vector of the grouping of SNPs for group penalties; SNPs within K1[i,i+1] are grouped together.
//'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
//'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
//'@param K0 a integer scalar for the number of SNPs/covariates that is not penalized.
//'@param tau a numeric scalar represents the element-wise soft-thresholding parameter for generating the LD.
//'@param eta a numeric scalar of the weights for incorporating prior information.
//'@param alpha a numeric scalar of the ridge penalty weights; alpha=1 corresponds to no ridge penalty.
//'@param eps a numeric scalar for the convergence criteria.
//'@param max_iter an integer scalar for the maximum iterations before the algorithm stops.
// [[Rcpp::export]]
List MaxLambda(Eigen::MatrixXd XtY, Eigen::MatrixXd tilde_beta, Eigen::SparseMatrix<double>& Sig,
                 Eigen::MatrixXd K1, Eigen::MatrixXd m, Eigen::MatrixXi& blk,int K0, double tau,
                 double eta, double alpha, double eps,int max_iter){
  // Lengths/dimensions
  int G = K1.rows() - 1;
  int P = Sig.cols();
  
  //Outcome
  Eigen::MatrixXd Lamb(G,1);
  Lamb.setZero(G,1);
  
  // Intermediate quantities
  Eigen::MatrixXd Xtbt(P,1);
  Xtbt.setZero(P,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd b(P,1);
  b.setZero(P,1);
  int tot_iter = 0;
  double shift=0, maxChange=0, len=0;
  
  // Initialization
  Xtbt=Sig*tilde_beta;
  
  while (tot_iter < max_iter) {
      
    tot_iter++;
    maxChange=0;
    
    // Update unpenalized covariates
    for (int j=0; j<K0; j++) {
      shift = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b)(0,0);
      if (fabs(shift) > maxChange) maxChange = fabs(shift);
      b(j,0) = shift + a(j,0);
    }
    
    // Check convergence
    a=b;
    if (maxChange <= eps) break;
  }
  
  for (int g=0; g<G; g++){
    
    int Kg = K1(g+1,0) - K1(g,0);
    // Calculate z
    Eigen::MatrixXd z(Kg,1);
    z.setZero(Kg,1);
    
    for (int j=K1(g,0); j<K1(g+1,0); j++) {
      z(j-K1(g,0),0) = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b)(0,0)+b(j,0);
    }
    Lamb(g,0) = Norm(z, Kg)*(1+eta)/alpha/m(g,0);
  }
  List result;
  //result["Sig"]=Sig;
  result["lambda.max"]=Lamb.maxCoeff();
  return(result);
}

//' BRIGHT estimation procedure
//'
//' Gradient-descent algorithm for optimizing the BRIGHT estimation procedure
//'
//'@param XtY a numeric verticle vector of the marginal correlation between genotype and outcome, can be recovered from GWAS summary statistics by function p2cor().
//'@param tilde_beta a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
//'@param X a numeric matrix of reference genotypes, default is from 1000 genome project.
//'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
//'@param K1 a numeric verticle vector of the grouping of SNPs for group penalties; SNPs within K1[i,i+1] are grouped together.
//'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
//'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
//'@param K0 a integer scalar for the number of SNPs/covariates that is not penalized.
//'@param penalty a integer scalar for chosing the penalties; 1 corresponds to LASSO, 2 corresponds to MCP, 3 corresponds to SCAD.
//'@param tau a numeric scalar represents the element-wise soft-thresholding parameter for generating the LD.
//'@param eta a numeric scalar of the weights for incorporating prior information.
//'@param alpha a numeric scalar of the ridge penalty weights; alpha=1 corresponds to no ridge penalty.
//'@param gamma a numeric scalar of the third parameter in firm and SCAD thresholding functions.
//'@param eps a numeric scalar for the convergence criteria.
//'@param max_iter an integer scalar for the maximum iterations before the algorithm stops.
//'@param dfmax an integer scalar for the maximum number of SNPs been selected before the algorithm stops.
//'@param gmax an integer scalar for the maximum number of groups been selected before the algorithm stops.
//'@param user a logical scalar for indicating whether lambda is user specified; if user=TRUE, then the iteration will start from the largest lambda, otherwise the iteration will start from the second largest lambda.
//'@return Write a list of results; beta, the coefficient estimate from BRIGHT; 
//'@return iter, the number of total iterations needed for the model to converge with each lambda; 
//'@return df total degree of freedom of the converged model with each lambda;
//'@return dev, the approximated deviance associated with each lambda.
// [[Rcpp::export]]
List gdfit_gaussian(Eigen::MatrixXd XtY, Eigen::MatrixXd tilde_beta, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXd lambda, 
                    Eigen::MatrixXd K1, Eigen::MatrixXd m, Eigen::MatrixXi& blk, int K0, int penalty, double tau,
                    double eta, double alpha, double gamma, double eps,int max_iter, int dfmax, int gmax, bool user) {
  
  // Lengths/dimensions
  int L = lambda.rows();
  int G = K1.rows() - 1;
  int P = Sig.cols();
  
  //Outcome
  List result;
  Eigen::MatrixXd b(P,L);
  b.setZero(P,L);
  Eigen::MatrixXd iter(L,1);
  iter.setZero(L,1);
  Eigen::MatrixXd df(L,1);
  df.setZero(L,1);
  Eigen::MatrixXd dev(L,1);
  dev.setZero(L,1);
  Eigen::MatrixXd dev2(L,1);
  dev2.setZero(L,1);
 
  // Intermediate quantities
  Eigen::MatrixXd Xtbt(P,1);
  Xtbt.setZero(P,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd e(G,1);
  e.setZero(G,1);
  int lstart=0, ng=0, nv=0, violations=0, tot_iter = 0;
  double shift=0, l1=0, l2=0, maxChange=0, maxChange2=0;
  
  // Initialization
  Xtbt=Sig*tilde_beta;
  
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
  }
  
  for (int l=lstart; l<L; l++) {
    //Rcout << l;
    //Rcout << "l\n";
    R_CheckUserInterrupt();
    if (l != 0) {
      a=b.block(0,l-1,P,1);
      
      // Check dfmax, gmax
      // predecide how many groups or parameters we want to select
      ng = 0;
      nv = 0;
      for (int g=0; g<G; g++) {
        if (a(K1(g,0),0) != 0) {
          ng++;
          nv = nv + (K1(g+1,0)-K1(g,0));
        }
      }
      if (ng > gmax || nv > dfmax || tot_iter == max_iter) {
        for (int ll=l; ll<L; ll++) iter(ll,0) = NA_INTEGER;
        Rcout<<"\nbreak through gmax and dfmax\n";
        Rcout<<ng;
        Rcout<<"\n";
        Rcout<<nv;
        Rcout<<"\n";
        Rcout<<gmax;
        Rcout<<"\n";
        Rcout<<dfmax;
        Rcout<<"\n";
        Rcout<<tot_iter;
        Rcout<<"\n";
        Rcout<<max_iter;
        break;
      }
    }
    
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        iter(l,0)++;
        tot_iter++;
        df(l,0)=0;
        maxChange=0;
        
        // Update unpenalized covariates
        for (int j=0; j<K0; j++) {
          shift = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b.block(0,l,P,1))(0,0);
          if (fabs(shift) > maxChange) maxChange = fabs(shift);
          b(j,l) = shift + a(j,0);
          df(l,0) += 1;
        }
        
        // Update penalized groups
        for (int g=0; g<G; g++) {
          if(e(g,0)==1){
            l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
            l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
            gd_gaussian(a, b, Sig, K1, Xtbt, XtY, df, l, P, g, penalty, l1, l2, gamma, eta, maxChange);
          }
        }
        
        // Check convergence
        a=b.block(0,l,P,1);
        if(maxChange==maxChange2) break;
        maxChange2=maxChange;
        //dev(l,0)=(1+eta)/2*(a.transpose()*Sig*a)(0,0)-(a.transpose()*(XtY+eta*Xtbt))(0,0);
        if (maxChange <= eps) break;
        
      }
      
      // Scan for violations
      violations = 0;
      
      for (int g=0; g<G; g++){
        if (e(g,0)==0) {
          l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
          l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
          gd_gaussian(a, b, Sig, K1, Xtbt, XtY, df, l, P, g, penalty, l1, l2, gamma, eta, maxChange);
          
          if (b(K1(g,0),l) != 0) {
	    //Rcout << g;
	    //Rcout << "added to active set\n";
            e(g,0) = 1;
            violations++;
          }
        }
      }
      
      if (violations==0){
        a=b.block(0,l,P,1);
        dev(l,0)=(1+eta)/2*(a.transpose()*Sig*a)(0,0)-(a.transpose()*(XtY+eta*Xtbt))(0,0);
	dev2(l,0)=(a.transpose()*Sig*a)(0,0)/2-(a.transpose()*(XtY))(0,0);
	//Rcout << (1+eta)/2*(a.transpose()*Sig*a)(0,0);
	//Rcout << "\n";
        //Rcout<<"active set done\n";
        //Rcout<<dev(l,0);
        //Rcout<<"Done\n";
        break;
      }
      a=b.block(0,l,P,1);
      //dev(l,0)=(1+eta)/2*(a.transpose()*Sig*a)(0,0)-(a.transpose()*(XtY+eta*Xtbt))(0,0);
      //Rcout<<"active set done\n";
    }
    
  }
  result["Beta"]=b;
  result["iter"]=iter;
  result["df"]=df;
  //result["Sig"]=Sig;
  result["dev"]=dev;
  result["dev2"]=dev2;
  return(result);
}


//Individual level fitting

// [[Rcpp::export]]
void gd_gaussiani(Eigen::MatrixXd & a, Eigen::MatrixXd & b, Eigen::MatrixXd & K1,
                  Eigen::MatrixXd & X, Eigen::MatrixXd & Yt, Eigen::MatrixXd & df,
                  int & l, int & P, int & g, int & penalty, double & lam1, double & lam2,
                  double & gamma, double & eta, double & maxChange) {
  
  int N = X.rows();
  int Kg = K1(g+1,0) - K1(g,0);
  // Calculate z
  Eigen::MatrixXd z(Kg,1);
  z.setZero(Kg,1);
  //for (int j=K1(g,0); j<K1(g+1,0); j++) {
  //  z(j-K1(g,0),0) = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b.block(0,l,P,1))(0,0)+b(j,l);
  
  /*
   double tmp=(Sig.col(K1(g,0)).transpose()*b.block(0,l,P,1))(0,0)+b(K1(g,0),l);
   if (NumericVector::is_na(tmp)){
   Rcout << j;
   Rcout << "j\n";
   Rcout << "first part\n";
   Rcout << Sig.col(j).transpose();
   Rcout << "\nfirst part end\n";
   Rcout << "Second part\n";
   Rcout << b.block(0,l,P,1).transpose();
   Rcout << "\nSecond part end\n";
   Rcout << "Third part\n";
   Rcout << b(j,l);
   Rcout << "\nThird part end\n";
   }
   */
  
  //}
  z = (X.block(0,K1(g,0),N,Kg).transpose()*(Yt-X*b.block(0,l,P,1))+b.block(K1(g,0),l,Kg,1))/N;
  double z_norm = Norm(z, Kg);
  
  if (NumericVector::is_na(z_norm)){
    Rcout << K1(g,0);
    Rcout << "K1(g,0)\n";
  }
  
  
  // Update b
  double len;
  if (penalty==1) len = S(z_norm, lam1) / (1+lam2);
  if (penalty==2) len = Ff(z_norm, lam1, lam2, gamma);
  if (penalty==3) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 || a(K1(g,0),0) != 0) {
    // If necessary, update beta and r
    for (int j=K1(g,0); j<K1(g+1,0); j++) {
      b(j,l) = len * z(j-K1(g,0),0) / z_norm;
      double shift = b(j,l)-a(j,0);
      /*Rcout << "shift\n";
       Rcout << shift;
       Rcout << "shift_end\n";
       if(fabs(shift)>10){
       Rcout << "len\n";
       Rcout << len;
       Rcout << "len_end\n";
       Rcout << "z\n";
       Rcout << z(j-K1(g,0),0);
       Rcout << "z_end\n";
       Rcout << "z_norm\n";
       Rcout << z_norm;
       Rcout << "z_norm_end\n";
       }*/
      
      if (fabs(shift) > maxChange) maxChange = fabs(shift);
    }
  }
  // Update df
  if (len > 0) df(l,0) += Kg * len / z_norm;
  
}

/*
 // [[Rcpp::export]]
 Eigen::SparseMatrix<double> Sigma_test(Eigen::MatrixXd& X, double tau, int P, int N){
 Eigen::SparseMatrix<double> Sig(P,P);
 Sig.reserve(Eigen::VectorXi::Constant(P,1000));
 Sigma(X, Sig, tau, P, N);
 return(Sig.col(1));
 }
 */

//' Maximum Lambda selection
//'
//' Calculate the maximum lambda, which corresponding to penalizing all coefficients to be zero
//'
//'@param XtY a numeric verticle vector of the marginal correlation between genotype and outcome, can be recovered from GWAS summary statistics by function p2cor().
//'@param tilde_beta a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
//'@param X a numeric matrix of reference genotypes, default is from 1000 genome project.
//'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
//'@param K1 a numeric verticle vector of the grouping of SNPs for group penalties; SNPs within K1[i,i+1] are grouped together.
//'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
//'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
//'@param K0 a integer scalar for the number of SNPs/covariates that is not penalized.
//'@param tau a numeric scalar represents the element-wise soft-thresholding parameter for generating the LD.
//'@param eta a numeric scalar of the weights for incorporating prior information.
//'@param alpha a numeric scalar of the ridge penalty weights; alpha=1 corresponds to no ridge penalty.
//'@param eps a numeric scalar for the convergence criteria.
//'@param max_iter an integer scalar for the maximum iterations before the algorithm stops.
// [[Rcpp::export]]
List MaxLambdai(Eigen::MatrixXd Y, Eigen::MatrixXd tilde_beta, Eigen::MatrixXd &X,
                Eigen::MatrixXd K1, Eigen::MatrixXd m,int K0, double tau,
                double eta, double alpha, double eps,int max_iter){
  // Lengths/dimensions
  int G = K1.rows() - 1;
  int P = X.cols();
  int N = X.rows();
  
  //Outcome
  Eigen::MatrixXd Lamb(G,1);
  Lamb.setZero(G,1);
  
  // Intermediate quantities
  Eigen::MatrixXd Xtbt(N,1);
  Xtbt.setZero(N,1);
  Eigen::MatrixXd Yt(N,1);
  Yt.setZero(N,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd b(P,1);
  b.setZero(P,1);
  int tot_iter = 0;
  double shift=0, maxChange=0, len=0;
  
  // Initialization
  Xtbt=X*tilde_beta;
  Yt=(Y+eta*Xtbt)/(1+eta);
  
  while (tot_iter < max_iter) {
    
    tot_iter++;
    maxChange=0;
    
    // Update unpenalized covariates
    for (int j=0; j<K0; j++) {
      shift = (X.block(0,j,N,1).transpose()*(Yt-(X*b.block(0,0,P,1))))(0,0)/N;
      if (fabs(shift) > maxChange) maxChange = fabs(shift);
      b(j,0) = shift + a(j,0);
    }
    
    // Check convergence
    a=b;
    if (maxChange <= eps) break;
  }
  
  for (int g=0; g<G; g++){
    
    int Kg = K1(g+1,0) - K1(g,0);
    // Calculate z
    Eigen::MatrixXd z(Kg,1);
    z.setZero(Kg,1);
    
    z = (X.block(0,K1(g,0),N,Kg).transpose()*(Yt-X*b.block(0,0,P,1))+b.block(K1(g,0),0,Kg,1))/N;
    
    Lamb(g,0) = Norm(z, Kg)*(1+eta)/alpha/m(g,0);
  }
  List result;
  //result["Sig"]=Sig;
  result["lambda.max"]=Lamb.maxCoeff();
  return(result);
}

//' BRIGHT estimation procedure
//'
//' Gradient-descent algorithm for optimizing the BRIGHT estimation procedure
//'
//'@param XtY a numeric verticle vector of the marginal correlation between genotype and outcome, can be recovered from GWAS summary statistics by function p2cor().
//'@param tilde_beta a numeric verticle vector of the prior information; in the case where prior information is a subspace of XtY, set the complement set to be zero.
//'@param X a numeric matrix of reference genotypes, default is from 1000 genome project.
//'@param lambda a numeric verticle vector of the LASSO penalty weights; it must be sorted with a decreasing order.
//'@param K1 a numeric verticle vector of the grouping of SNPs for group penalties; SNPs within K1[i,i+1] are grouped together.
//'@param m a numeric verticle vector of the multiplicative factor to adjust for the number of SNPs in each group, the default is the square root of number of elements in each group.
//'@param blk a numeric verticle vector indicating the grouping of SNPs in each block, can be generated through LD() function.
//'@param K0 a integer scalar for the number of SNPs/covariates that is not penalized.
//'@param penalty a integer scalar for chosing the penalties; 1 corresponds to LASSO, 2 corresponds to MCP, 3 corresponds to SCAD.
//'@param tau a numeric scalar represents the element-wise soft-thresholding parameter for generating the LD.
//'@param eta a numeric scalar of the weights for incorporating prior information.
//'@param alpha a numeric scalar of the ridge penalty weights; alpha=1 corresponds to no ridge penalty.
//'@param gamma a numeric scalar of the third parameter in firm and SCAD thresholding functions.
//'@param eps a numeric scalar for the convergence criteria.
//'@param max_iter an integer scalar for the maximum iterations before the algorithm stops.
//'@param dfmax an integer scalar for the maximum number of SNPs been selected before the algorithm stops.
//'@param gmax an integer scalar for the maximum number of groups been selected before the algorithm stops.
//'@param user a logical scalar for indicating whether lambda is user specified; if user=TRUE, then the iteration will start from the largest lambda, otherwise the iteration will start from the second largest lambda.
//'@return Write a list of results; beta, the coefficient estimate from BRIGHT; 
//'@return iter, the number of total iterations needed for the model to converge with each lambda; 
//'@return df total degree of freedom of the converged model with each lambda;
//'@return dev, the approximated deviance associated with each lambda.
// [[Rcpp::export]]
List gdfit_gaussiani(Eigen::MatrixXd Y, Eigen::MatrixXd &X, Eigen::MatrixXd tilde_beta, Eigen::MatrixXd lambda, 
                     Eigen::MatrixXd K1, Eigen::MatrixXd m, int K0, int penalty, double tau,
                     double eta, double alpha, double gamma, double eps,int max_iter, int dfmax, int gmax, bool user) {
  
  // Lengths/dimensions
  int L = lambda.rows();
  int G = K1.rows() - 1;
  int P = X.cols();
  int N = X.rows();
  
  //Outcome
  List result;
  Eigen::MatrixXd b(P,L);
  b.setZero(P,L);
  Eigen::MatrixXd iter(L,1);
  iter.setZero(L,1);
  Eigen::MatrixXd df(L,1);
  df.setZero(L,1);
  Eigen::MatrixXd dev(L,1);
  dev.setZero(L,1);
  Eigen::MatrixXd dev2(L,1);
  dev2.setZero(L,1);
  
  // Intermediate quantities
  Eigen::MatrixXd Xtbt(N,1);
  Xtbt.setZero(N,1);
  Eigen::MatrixXd Yt(N,1);
  Yt.setZero(N,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd e(G,1);
  e.setZero(G,1);
  int lstart=0, ng=0, nv=0, violations=0, tot_iter = 0;
  double shift=0, l1=0, l2=0, maxChange=0, maxChange2=0;
  
  // Initialization
  Xtbt=X*tilde_beta;
  Yt=(Y+eta*Xtbt)/(1+eta);
  
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
  }
  
  for (int l=lstart; l<L; l++) {
    Rcout << l;
    Rcout << "l\n";
    R_CheckUserInterrupt();
    if (l != 0) {
      a=b.block(0,l-1,P,1);
      
      // Check dfmax, gmax
      // predecide how many groups or parameters we want to select
      ng = 0;
      nv = 0;
      for (int g=0; g<G; g++) {
        if (a(K1(g,0),0) != 0) {
          ng++;
          nv = nv + (K1(g+1,0)-K1(g,0));
        }
      }
      if (ng > gmax || nv > dfmax || tot_iter == max_iter) {
        for (int ll=l; ll<L; ll++) iter(ll,0) = NA_INTEGER;
        Rcout<<"\nbreak through gmax and dfmax\n";
        Rcout<<ng;
        Rcout<<"\n";
        Rcout<<nv;
        Rcout<<"\n";
        Rcout<<gmax;
        Rcout<<"\n";
        Rcout<<dfmax;
        Rcout<<"\n";
        Rcout<<tot_iter;
        Rcout<<"\n";
        Rcout<<max_iter;
        break;
      }
    }
    
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        iter(l,0)++;
        tot_iter++;
        df(l,0)=0;
        maxChange=0;
        
        // Update unpenalized covariates
        for (int j=0; j<K0; j++) {
          shift = (X.block(0,j,N,1).transpose()*(Yt-(X*b.block(0,l,P,1))))(0,0)/N;
          if (fabs(shift) > maxChange) maxChange = fabs(shift);
          b(j,l) = shift + a(j,0);
          df(l,0) += 1;
        }
        
        // Update penalized groups
        for (int g=0; g<G; g++) {
          if(e(g,0)==1){
            l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
            l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
            gd_gaussiani(a, b, K1, X, Yt, df, l, P, g, penalty, l1, l2, gamma, eta, maxChange);
          }
        }
        
        // Check convergence
        a=b.block(0,l,P,1);
        if(maxChange==maxChange2) break;
        maxChange2=maxChange;
        //dev(l,0)=(1+eta)/2*(a.transpose()*Sig*a)(0,0)-(a.transpose()*(XtY+eta*Xtbt))(0,0);
        if (maxChange <= eps) break;
        
      }
      
      // Scan for violations
      violations = 0;
      
      for (int g=0; g<G; g++){
        if (e(g,0)==0) {
          l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
          l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
          gd_gaussiani(a, b, K1, X, Yt, df, l, P, g, penalty, l1, l2, gamma, eta, maxChange);
          
          if (b(K1(g,0),l) != 0) {
            //Rcout << g;
            //Rcout << "added to active set\n";
            e(g,0) = 1;
            violations++;
          }
        }
      }
      
      if (violations==0){
        a=b.block(0,l,P,1);
        dev(l,0)=((Yt-X*a).transpose()*(Yt-X*a))(0,0)/N;
        dev2(l,0)=((Y-X*a).transpose()*(Y-X*a))(0,0)/N;
        //Rcout << (1+eta)/2*(a.transpose()*Sig*a)(0,0);
        //Rcout << "\n";
        //Rcout<<"active set done\n";
        //Rcout<<dev(l,0);
        //Rcout<<"Done\n";
        break;
      }
      a=b.block(0,l,P,1);
      //dev(l,0)=(1+eta)/2*(a.transpose()*Sig*a)(0,0)-(a.transpose()*(XtY+eta*Xtbt))(0,0);
      //Rcout<<"active set done\n";
    }
    
  }
  result["Beta"]=b;
  result["iter"]=iter;
  result["df"]=df;
  //result["Sig"]=Sig;
  result["dev"]=dev;
  result["dev2"]=dev2;
  return(result);
}

