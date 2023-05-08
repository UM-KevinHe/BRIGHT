#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>

// Cross product of the jth column of x with y
double crossprod(double *x, double *y, int n, int j) {
  double val = 0;
  int nn = n*j;
  for (int i=0; i<n; i++) val += x[nn+i] * y[i];
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of x
double sum(double *x, int n) {
  double val = 0;
  for (int i=0; i<n; i++) val += x[i];
  return(val);
}

// Max of x
double max(double *x, int n) {
  double val = x[0];
  for (int i=1; i<n; i++) {
    if (x[i] > val) val = x[i];
  }
  return(val);
}

// Gaussian loss
double gLoss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Pr(y=1) for binomial
double p_binomial(double eta) {
  if (eta > 10) {
    return(1);
  } else if (eta < -10) {
    return(0);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}

// Euclidean norm
double norm(double *x, int p) {
  double x_norm = 0;
  for (int j=0; j<p; j++) x_norm = x_norm + pow(x[j],2);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

// Soft-thresholding operator
double S(double z, double l) {
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

// Firm-thresholding operator
double F(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(1+l2-1/gamma));
  else return(z/(1+l2));
}

// SCAD-modified firm-thresholding operator
double Fs(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(1+l2));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(1-1/(gamma-1)+l2));
  else return(z/(1+l2));
}

// MCP penalty
double MCP(double theta, double l, double a) {
  theta = fabs(theta);
  if (theta <= a*l) return(l*theta - pow(theta,2)/(2*a));
  else return(a*pow(l,2)/2);
}

// MCP penalization rate
double dMCP(double theta, double l, double a) {
  theta = fabs(theta);
  if (theta < a*l) return(l-theta/a);
  else return(0);
}


SEXP standardize(SEXP X_) {
  // Declarations
  int n = nrows(X_);
  int p = ncols(X_);
  SEXP XX_, c_, s_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  for (int j=0; j<(n*p); j++) REAL(XX_)[j] = 0;
  PROTECT(c_ = allocVector(REALSXP, p));
  for (int j=0; j<p; j++) REAL(c_)[j] = 0;
  PROTECT(s_ = allocVector(REALSXP, p));
  for (int j=0; j<p; j++) REAL(s_)[j] = 0;
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *c = REAL(c_);
  double *s = REAL(s_);
  
  for (int j=0; j<p; j++) {
    // Center
    c[j] = 0;
    for (int i=0; i<n; i++) {
      c[j] += X[j*n+i];
    }
    c[j] = c[j] / n;
    for (int i=0; i<n; i++) XX[j*n+i] = X[j*n+i] - c[j];
    
    // Scale
    s[j] = 0;
    for (int i=0; i<n; i++) {
      s[j] += pow(XX[j*n+i], 2);
    }
    s[j] = sqrt(s[j]/n);
    for (int i=0; i<n; i++) XX[j*n+i] = XX[j*n+i]/s[j];
  }
  
  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_);
  SET_VECTOR_ELT(res, 1, c_);
  SET_VECTOR_ELT(res, 2, s_);
  UNPROTECT(4);
  return(res);
}


SEXP maxprod(SEXP X_, SEXP y_, SEXP K_, SEXP m_) {
  
  // Declarations
  int n = nrows(X_);
  int J = length(K_)-1;
  SEXP zmax;
  PROTECT(zmax = allocVector(REALSXP, 1));
  REAL(zmax)[0] = 0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *K = INTEGER(K_);
  
  for (int g=0; g<J; g++) {
    for (int j=K[g]; j<K[g+1]; j++) {
      double z = fabs(crossprod(X, y, n, j) / m[g]);
      if (z > REAL(zmax)[0]) REAL(zmax)[0] = z;
    }
  }
  
  // Return
  UNPROTECT(1);
  return(zmax);
}


SEXP maxgrad(SEXP X_, SEXP y_, SEXP K_, SEXP m_) {
  
  // Declarations
  int n = nrows(X_);
  int J = length(K_)-1;
  SEXP zmax;
  PROTECT(zmax = allocVector(REALSXP, 1));
  REAL(zmax)[0] = 0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *K = INTEGER(K_);
  
  for (int g=0; g<J; g++) {
    int Kg = K[g+1]-K[g];
    double *Z = Calloc(Kg, double);
    for (int j=K[g]; j<K[g+1]; j++) {
      Z[j-K[g]] = crossprod(X, y, n, j);
    }
    double z = norm(Z, Kg) / m[g];
    if (z > REAL(zmax)[0]) REAL(zmax)[0] = z;
    Free(Z);
  }
  
  // Return
  UNPROTECT(1);
  return(zmax);
}


void gd_glm(double *b, double *x, double *r, double v, double *eta, int g, int *K1, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a, double *maxChange) {
  
  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + a[j]; // Discrete this r will be t*n dimensional
  double z_norm = norm(z,K);
  
  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(v * z_norm, lam1) / (v * (1 + lam2));// why are we multiplying v? MM doesn't need to argue 2 derivative
  if (strcmp(penalty, "grMCP")==0) len = F(v * z_norm, lam1, lam2, gamma) / v;
  if (strcmp(penalty, "grSCAD")==0) len = Fs(v * z_norm, lam1, lam2, gamma) / v;
  if (len != 0 || a[K1[g]] != 0) {// a[K1[g]]!=0, len = 0, it will panelize g to be 0
    // If necessary, update b and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      b[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = b[l*p+j]-a[j];
      if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
      for (int i=0; i<n; i++) {
        double si = shift*x[j*n+i];
        r[i] -= si;
        eta[i] += si;
      }
    }
  }
  
  // Update df
  if (len > 0) REAL(df)[l] += K * len / z_norm;
  Free(z);
}



SEXP gdfit_glm(SEXP X_, SEXP y_, SEXP family_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_) {
  
  // Lengths/dimensions
  int n = length(y_);
  int L = length(lambda);
  int J = length(K1_) - 1;
  int p = length(X_)/n;
  
  // Pointers
  double *X = REAL(X_);
  double *y = REAL(y_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  const char *family = CHAR(STRING_ELT(family_, 0));
  int *K1 = INTEGER(K1_);
  int K0 = INTEGER(K0_)[0];
  double *lam = REAL(lambda);
  double alpha = REAL(alpha_)[0];
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int warn = INTEGER(warn_)[0];
  int user = INTEGER(user_)[0];
  
  // Outcome
  SEXP res, beta0, beta, iter, df, Dev;
  PROTECT(res = allocVector(VECSXP, 5));
  PROTECT(beta0 = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(beta0)[i] = 0;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(Dev = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(Dev)[i] = 0;
  double *b0 = REAL(beta0);
  double *b = REAL(beta);
  
  // Intermediate quantities
  double a0 = 0; // Beta0 from previous iteration
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *eta = Calloc(n, double);
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int lstart, ng, nv, violations;
  double shift, l1, l2, mu, v, maxChange;
  
  // Initialization
  double ybar = sum(y, n)/n;
  double nullDev = 0;
  if (strcmp(family, "binomial")==0) {
    a0 = b0[0] = log(ybar/(1-ybar));
    for (int i=0; i<n; i++) nullDev -= 2*y[i]*log(ybar) + 2*(1-y[i])*log(1-ybar);
  } else if (strcmp(family, "poisson")==0) {
    a0 = b0[0] = log(ybar);
    for (int i=0;i<n;i++) {
      if (y[i]!=0) nullDev += 2*(y[i]*log(y[i]/ybar) + ybar - y[i]);
      else nullDev += 2*ybar;
    }
  }
  for (int i=0; i<n; i++) eta[i] = a0;
  
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Dev)[0] = nullDev;
  }
  
  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      a0 = b0[l-1];
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];
      
      // Check dfmax, gmax
      // predecide how many groups or parameters we want to select
      ng = 0;
      nv = 0;
      for (int g=0; g<J; g++) {
        if (a[K1[g]] != 0) {
          ng++;
          nv = nv + (K1[g+1]-K1[g]);
        }
      }
      if (ng > gmax || nv > dfmax || tot_iter == max_iter) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }
    }
    
    
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) { //Two loop add up to max_iter, can be problematic, because its combined together, inner loop can go to far
        INTEGER(iter)[l]++;
        tot_iter++;
        
        // Approximate L
        REAL(Dev)[l] = 0;
        if (strcmp(family, "binomial")==0) {
          v = 0.25; //second derivative approximation, MM algorithm, second derivative is vector so the computing speed is not reduced very much
          for (int i=0; i<n; i++) {
            mu = p_binomial(eta[i]);
            r[i] = (y[i] - mu) / v;// L2^{-1}*L1 (newton raphson direction) not Score function, r is residual of pseudo outcome, z-xbeta^hat_{-j}
            REAL(Dev)[l] -= 2*y[i]*log(mu)+2*(1-y[i])*log(1-mu); //Should be changed for our project
          }
        } else if (strcmp(family, "poisson")==0) {
          v = exp(max(eta, n));
          for (int i=0; i<n; i++) {
            mu = exp(eta[i]);
            r[i] = (y[i] - mu)/v;
            if (y[i]!=0) REAL(Dev)[l] += 2*(y[i]*log(y[i]/mu) + mu - y[i]);
            else REAL(Dev)[l] += 2*mu;
          }
        }
        
        // Check for saturation
        if (REAL(Dev)[l]/nullDev < .01) {//early stopping? print(REAL(Dev)[l]/nullDev)
          if (warn) warning("Model saturated; exiting...");
          for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }
        
        // Update intercept
        // Coordinate descent algorithm, speed may be faster
        shift = sum(r, n)/n;// 1 step newton raphson?
        b0[l] = shift + a0;
        for (int i=0; i<n; i++) {
          r[i] -= shift;
          eta[i] += shift;
        }
        REAL(df)[l] = 1;
        maxChange = fabs(shift);
        
        // Update unpenalized covariates
        for (int j=0; j<K0; j++) {
          shift = crossprod(X, r, n, j)/n;
          if (fabs(shift) > maxChange) maxChange = fabs(shift);
          b[l*p+j] = shift + a[j];
          for (int i=0; i<n; i++) {
            double si = shift * X[n*j+i];
            r[i] -= si;// Pseudo outcome
            eta[i] += si;
          }
          REAL(df)[l]++;
        }
        
        // Update penalized groups
        for (int g=0; g<J; g++) {
          l1 = lam[l] * m[g] * alpha;
          l2 = lam[l] * m[g] * (1-alpha);
          if (e[g]) gd_glm(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, &maxChange);
        }
        
        // Check convergence
        a0 = b0[l];
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }
      
      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
        if (e[g]==0) {
          l1 = lam[l] * m[g] * alpha;
          l2 = lam[l] * m[g] * (1-alpha);
          gd_glm(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, &maxChange);
          if (b[l*p+K1[g]] != 0) {
            e[g] = 1;
            violations++;
          }
        }
      }
      
      if (violations==0) break;
      a0 = b0[l];
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  Free(a);
  Free(r);
  Free(e);
  Free(eta);
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, df);
  SET_VECTOR_ELT(res, 4, Dev);
  UNPROTECT(6);
  return(res);
}



