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
double F(double z, double l1, double l2, double gamma) {
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
int Sigma(Eigen::MatrixXd& X, Eigen::SparseMatrix<double>& Sig, double tau, int P, int N){
  for(int i=0;i<P;i++){
    Sig.insert(i,i) = 1.0;
    for(int j=i+1;j<P;j++){
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

// [[Rcpp::export]]
void gd_gaussian(Eigen::MatrixXd & a, Eigen::MatrixXd & b, Eigen::MatrixXd & Sig, Eigen::MatrixXd & K1,
                 Eigen::MatrixXd & Xtbt, Eigen::MatrixXd & XtY, Eigen::MatrixXd & df,
                 int & l, int & P, int & g, int & penalty, double & lam1, double & lam2,
                 double & gamma, double & eta, double & maxChange) {
  
  int Kg = K1(g+1,0) - K1(g,0);
  // Calculate z
  Eigen::MatrixXd z(Kg,0);
  z.setZero(Kg,0);
  for (int j=K1(g,0); j<K1(g+1,0); j++) z(j-K1(g,0),0) = (XtY(j,0)+eta*Xtbt(j,0))/(1+eta)-(Sig.col(j).transpose()*b.block(0,l,P,1))(0,0)+b(j,0);
  double z_norm = norm(z, Kg);
  
  // Update b
  double len;
  if (penalty==1) len = S(z_norm, lam1) / (1+lam2);
  if (penalty==2) len = F(z_norm, lam1, lam2, gamma);
  if (penalty==3) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 || a(K1(g,0),0) != 0) {
    // If necessary, update beta and r
    for (int j=K1(g,0); j<K1(g+1,0); j++) {
      b(j,l) = len * z(j-K1(g,0),0) / z_norm;
      double shift = b(j,l)-a(j,0);
      if (fabs(shift) > maxChange) maxChange = fabs(shift);
    }
  }
  // Update df
  if (len > 0) df(l,0) += Kg * len / z_norm;
  Free(z);
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
List gdfit_gaussian(Eigen::MatrixXd XtY, Eigen::MatrixXd tilde_beta, Eigen::MatrixXd &X, Eigen::MatrixXd lambda, double eta,Eigen::MatrixXd K1, int K0, 
                    double lam_max, double alpha, double eps, 
                    int max_iter, double group_multiplier, 
                    int dfmax, int gmax, bool user) {
  
  // Lengths/dimensions
  int L = lambda.rows();
  int G = K1.rows() - 1;
  int P = X.cols();
  int N = X.rows();
  
  //Outcome
  List result;
  Eigen::MatrixXd beta(P,L);
  beta.setZero(P,L);
  Eigen::MatrixXd iter(L,1);
  iter.setZero(L,1);
  Eigen::MatrixXd df(L,1);
  df.setZero(L,1);
  Eigen::MatrixXd Dev(L,1);
  Dev.setZero(L,1);
  Eigen::MatrixXd b = beta;
  
  // Intermediate quantities
  Eigen::SparseMatrix<double> Sig(P,P);
  Sig.reserve(Eigen::VectorXi::Constant(P,1000));
  Eigen::MatrixXd Xtbt(P,1);
  Xtbt.setZero(P,1);
  Eigen::MatrixXd a(P,1);
  a.setZero(P,1);
  Eigen::MatrixXd e(G,1);
  e.setZero(G,1);
  int lstart=0, ng=0, nv=0, violations=0, tot_iter = 0;
  double shift=0, l1=0, l2=0, v=0, si=0, maxChange=0, z_norm=0, len=0;
  
  // Initialization
  
  Xtr=XtY;
  Sigma(X, Sig, tau, P, N);
  Xtbt=Sig*tilde_beta;
  
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
  }
  
  for (int l=lstart; l<L; l++) {
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
          l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
          l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
          if(e(g,0)==1){
            gd_gaussian(a, b, Sig, K1, Xtbt, XtY, df, l, P, g, penalty, l1, l2, gamma=0, eta, maxChange);
          }
        }
        
        // Check convergence
        a=b.block(0,l,P,1);
        if (maxChange <= eps) break;
        
      }
      
      // Scan for violations
      violations = 0;
      
      for (int g=0; g<G; g++){
        if (e(g,0)==0) {
          l1 = lambda(l,0) * m(g,0) * alpha/(1+eta);
          l2 = lambda(l,0) * m(g,0) * (1-alpha)/(1+eta);
          gd_gaussian(a, b, Sig, K1, Xtbt, XtY, df, l, P, g, penalty, l1, l2, gamma=0, eta, maxChange);
          
          if (b(K1(g,0),l) != 0) {
            e(g,0) = 1;
            violations++;
          }
        }
      }
      
      if (violations==0) break;
      a=b.block(0,l,P,1);
    }
    
  }
  result["Betat"]=b0;
  result["Betav"]=b;
  result["iter"]=iter;
  result["df"]=df;
  result["Dev"]=Dev;
  result["eta"]=eta;
  result["residual"]=r;
  return(result);
}