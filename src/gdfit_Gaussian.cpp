#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <string.h>
using namespace Rcpp;

//' Second norm
//'
//' Calculate second norm of a vertical vector
//'
//'@param x a numeric column vector to calculate the second norm
//'@param p a integer represents the dimension of x
//'@return ||x||_2
// [[Rcpp::export]]
double norm(Eigen::MatrixXd x, int p) {
  double x_norm = 0;
  for (int j=0; j<p; j++) x_norm = x_norm + x(j,0)*x(j,0);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

//' Soft-thresholding
//'
//' Calculate the soft-thresholding result
//'
//'@param z a numeric scalar to be soft-thresholded
//'@param l a numeric scalar for the soft-thresholding parameter
//'@return S(z,l)
// [[Rcpp::export]]
double S(double z, double l) {
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

//' Firm-thresholding
//'
//' Calculate the Firm-thresholding result
//'
//'@param z a numeric scalar to be Firm-thresholded
//'@param l1 a numeric scalar for the first Firm-thresholding parameter
//'@param l2 a numeric scalar for the second Firm-thresholding parameter
//'@param gamma a numeric scalar for the third Firm-thresholding parameter
//'@return F(z,l1,l2)
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
//' Calculate the SCAD-thresholding result
//'
//'@param z a numeric scalar to be SCAD-thresholded
//'@param l1 a numeric scalar for the first SCAD-thresholding parameter
//'@param l2 a numeric scalar for the second SCAD-thresholding parameter
//'@param gamma a numeric scalar for the third SCAD-thresholding parameter
//'@return Fs(z,l1,l2)
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

// [[Rcpp::export]]
int Sigma_LD(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, Eigen::MatrixXi& blk, double tau, int P, int N){
  int lblk=blk.rows();
  
  
  typedef Eigen::Triplet<double> T;
  std::vector< T > TL;
  TL.reserve(10000000);
  
  for(int k=0;k<lblk-1;k++){
    Eigen::MatrixXd mat=(X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)).transpose()*X.block(0,blk(k,0),N,blk(k+1,0)-blk(k,0)))/N;
      for(int i=blk(k,0);i<blk(k+1,0);i++){
      Rcout << i;
      Rcout << "\n";
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
  return(0);
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

  double z_norm = norm(z, Kg);
  
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

// [[Rcpp::export]]
List MaxLambda(Eigen::MatrixXd XtY, Eigen::MatrixXd tilde_beta, Eigen::MatrixXd &X,
                 Eigen::MatrixXd K1, Eigen::MatrixXd m, Eigen::MatrixXi& blk,int K0, double tau,
                 double eta, double alpha, double eps,int max_iter){
  // Lengths/dimensions
  int G = K1.rows() - 1;
  int P = X.cols();
  int N = X.rows();
  
  //Outcome
  Eigen::MatrixXd Lamb(G,1);
  Lamb.setZero(G,1);
  
  // Intermediate quantities
  Eigen::SparseMatrix<double, 0> Sig(P,P);
  Sig.reserve(Eigen::VectorXi::Constant(P,3000));
  Eigen::MatrixXd Xtbt(P,1);
  Xtbt.setZero(P,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd b(P,1);
  b.setZero(P,1);
  int tot_iter = 0;
  double shift=0, maxChange=0, len=0;
  
  // Initialization
  Sigma_LD(X, Sig, blk, tau, P, N);
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
    Lamb(g,0) = norm(z, Kg)*(1+eta)/alpha/m(g,0);
  }
  List result;
  //result["Sig"]=Sig;
  result["lambda.max"]=Lamb.maxCoeff();
  return(result);
}

// [[Rcpp::export]]
List gdfit_gaussian(Eigen::MatrixXd XtY, Eigen::MatrixXd tilde_beta, Eigen::MatrixXd &X, Eigen::MatrixXd lambda, 
                    Eigen::MatrixXd K1, Eigen::MatrixXd m, Eigen::MatrixXi& blk, int K0, int penalty, double tau,
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
  Eigen::SparseMatrix<double, 0> Sig(P,P);
  Sig.reserve(Eigen::VectorXi::Constant(P,3000));
  Eigen::MatrixXd Xtbt(P,1);
  Xtbt.setZero(P,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd e(G,1);
  e.setZero(G,1);
  int lstart=0, ng=0, nv=0, violations=0, tot_iter = 0;
  double shift=0, l1=0, l2=0, maxChange=0, maxChange2=0;
  
  // Initialization

  Sigma_LD(X, Sig, blk, tau, P, N);
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
