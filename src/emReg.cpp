// [[Rcpp::depends(BH)]]
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include "bzinb.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
using namespace std;
using namespace Rcpp;
#define ITER_ALLOWANCE 100   // number of iterations allowed after finding the peak
// #define EPSILON1 1e-7  //for E-steps -> tol
#define EPSILON2 1e-7  //for M-steps
#define EPSILON3 1e-4  //for M-steps (gamma)
// #define DEBUG4
// #define DEBUG5
// #define DEBUG7
//#define DEBUG8
//#define DEBUG8g

// 3. EM
// [[Rcpp::export]]
List emReg(NumericVector& param2, IntegerVector &xvec, IntegerVector &yvec, 
           NumericMatrix& ZZ, NumericMatrix& WW,
           int &pZ, int &pW,
           int &n, int &se, int &maxiter, double &tol, int showFlag,
           IntegerVector &zi)
{
  // zi = (indicators for gamma1, gamma2, gamma3, zi.index (in 0:4), dim_pi)
  NumericVector param = clone(param2);
  int dim_pi = zi[4]; // dim_pi = 0 (for zi = 0), 1 (for zi = 1,2,3), 3 (for zi = 4)
  int dim_param = 3 + 2 * pZ + dim_pi * pW;
  int dim_gam = dim_pi > 0 ? (pW * dim_pi): 1;
  
  NumericVector alpha(3);
  NumericVector eta1(pZ);
  NumericVector eta2(pZ);
  NumericVector gamma1(zi[0]? pW: 1, 0.0);
  NumericVector gamma2(zi[1]? pW: 1, 0.0);
  NumericVector gamma3(zi[2]? pW: 1, 0.0);
  double error;

  double param_diff = 1.0;
  NumericVector param_old(dim_param);
  //double idgam[3];
  //double lb[1];  // replaced by eta1
  NumericMatrix expt(12, n);
  NumericVector exptBar(12);
  NumericVector b1(n);
  NumericVector b2(n);
  NumericVector p1(n, 1.0);
  NumericVector p2(n, 0.0);
  NumericVector p3(n, 0.0);
  NumericVector p4(n, 0.0);
  NumericVector b21(n); // b2 / b1
  NumericVector Zbar(pZ);
  NumericVector epsilon(pZ);
  NumericVector gamma(dim_gam, 0.0);
  NumericVector gammaNew(dim_gam, 0.0);
  NumericVector AEnew(pZ + 3);
  NumericVector epsNew(pZ);
  IntegerVector freq(n, 1L); // For compatibility use only

  // alpha & eta1
  arma::mat V_eta1(pZ + 3, pZ + 3);
  arma::mat score(pZ + 3, 1);
  
  // epsilon
  arma::mat V_eps(pZ, pZ);
  arma::mat score_eps(pZ, 1);
  
  // gamma
  arma::mat V_gam(dim_gam, dim_gam);
  arma::mat score_gam(dim_gam, 1);
  
  // NumericVector expt(12, 0.0);
  NumericVector s_i_abp(8, 0.0);
  NumericVector s_i(dim_param, 0.0);
  NumericVector info(se ? dim_param * dim_param : 1L, 0.0);  // if se = 1, 8 x 8 matrix, o/w a scalar zero.
  IntegerVector iter(1, 1L);
  IntegerVector nonconv(1, 0L);
  NumericVector trajectory(maxiter + 1L, 0.0); 
  NumericVector expt_i(12, 0.0);

  // storage for max_likelihood.
  double param_max[dim_param];
  int iter_max = 0;
  arma::mat expt_max(12, n);
  double ll_max = 0.0;
  
  for (int i = 0; i < dim_param; i++) {
    param[i] = (double) param2[i];  //initializing param with param2
  }
  // assigning each parameter
  int offset = 0;
  for (int i = 0; i < 3; i++) alpha[i] = (double) param[offset + i];
  offset = 3;
  for (int i = 0; i < pZ; i++) eta1[i] = (double) param[offset + i];
  offset += pZ;
  for (int i = 0; i < pZ; i++) eta2[i] = (double) param[offset + i];
  offset += pZ;
  if (zi[0]) {
    for (int i = 0; i < pW; i++) gamma1[i] = (double) param[offset + i];
    offset += pW;
  }
  if (zi[1]) {
    for (int i = 0; i < pW; i++) gamma2[i] = (double) param[offset + i];
    offset += pW;
  } 
  if (zi[2]) {
    for (int i = 0; i < pW; i++) gamma3[i] = (double) param[offset + i];
  }
  
  for (int j = 0; j < pZ; j++) {
    Zbar[j] = 0.0;
    for (int i = 0; i < n; i++) {
      Zbar[j] += ZZ(i, j);
    }
    Zbar[j] /= n;
  }
  

  //cout << "maxiter = " << maxiter << " iter = " << iter[0] << endl;  
  while(maxiter >= iter[0] && param_diff > tol)
  {
// if (iter[0]> 3528) Rcout << "iter = " << iter[0] << endl;
    if (iter[0] < 3) {
      ll_max = exptBar[0];
    }
    // for (int i = 0;i < dim_param;i++)
    // {
    //   param_old[i] = param[i];
    //   cout << "param [" << i << "] = " << param[i] << " ";
    // }


    if (iter[0] <= 1L) {
      // reparametrization
      logLin(ZZ, eta1, n, pZ, 1, b1);
      logLin(ZZ, eta2, n, pZ, 1, b2);
      for (int i = 0; i < n; i++) b21[i] = b2[i] / b1[i];

      if (zi[0]) logLin(WW, gamma1, n, pW, 1, p2);
      if (zi[1]) logLin(WW, gamma2, n, pW, 1, p3);
      if (zi[2]) logLin(WW, gamma3, n, pW, 1, p4);
      for (int i= 0; i < n; i++) { // normalization and p4
        p1[i] = 1.0;
        if (zi[0]) p1[i] += p2[i];
        if (zi[1]) p1[i] += p3[i];
        if (zi[2]) p1[i] += p4[i];
        
        if (zi[0]) p2[i] /= p1[i];
        if (zi[1]) p3[i] /= p1[i];
        if (zi[2]) p4[i] /= p1[i];
        p1[i] = 1.0 / p1[i];
      }
    }
    
    /* Rcout << "b1 = ";
    for (int i=0; i<n; i++) {
      Rcout << b1[i] << " ";  
    }
    Rcout << endl << "b2 = "; 
    for (int i=0; i<n; i++) {
      Rcout << b2[i] << " ";  
    }
    Rcout << endl;
    Rcout << "pi1 = ";
    for (int i=0; i<n; i++) {
      Rcout << p1[i] << " ";  
    }
    Rcout << endl; */
    // Rcout << "E-step" << endl;
    dBvZINB_Expt_mat(xvec, yvec, ZZ, WW, n, pZ, pW,
                     alpha, b1, b2, p1, p2, p3, p4, 
                     expt, s_i, s_i_abp, info, 0, zi[3], expt_i);

// Rcout << "expt_mat = " << endl;
// for (int i = 0; i < n; i++) {
//   if (xvec[i]==0 & yvec[i] > 0) {
//     Rcout << i << ": ";
//     for (int j = 0; j < 12; j++) {
//       Rcout << expt(j, i) << " ";
//     }
//     Rcout << endl;
//   }
// }
// Rcout << endl;
    
    // exptBar.
    for (int i = 0; i < 12; i++) {
      exptBar[i] = 0.0;
      for (int j = 0; j < n; j++) {
        exptBar[i] += expt(i, j);
      }
    }
    for (int i = 1; i < 12; i++) {
      exptBar[i] /= n; // exclude normalizing the likelihood
    }
  
// Rcout << "exptBar" <<endl;
// for (int i = 0; i < 12; i++) {
//   Rcout << exptBar[i] << " ";
// }
// Rcout << endl << "expt_i" <<endl;
// for (int j = 0; j < n; j++) {
//   for (int i = 0; i < 12; i++) {
//     Rcout << expt(i, j) << " ";
//   }
//   Rcout << endl;
// }
    
    
    if ((showFlag > 0) & (iter[0] >= showFlag)) {
      Rcout << "iter " << iter[0]  << " lik " << exptBar[0] <<
        ", a(0,1,2) = (" << param[0] << ", " << param[1] << ", " << param[2] << 
          "), eta1 = (";
      for (int i = 0; i < pZ; i ++ ) 
        Rcout << param[3 + i] << ", ";
      Rcout << "), eta2 = (";
      for (int i = 0; i < pZ; i ++ ) 
        Rcout << param[3 + pZ + i] << " ";
      offset = 3 + 2 * pZ;
      if (zi[0]) {
        Rcout << "), gamma1 = (";
        for (int i = 0; i < pW; i ++ ) 
          Rcout << param[offset + i] << " ";
        offset += pW;
      }
      if (zi[1]) {
        Rcout << "), gamma2 = (";
        for (int i = 0; i < pW; i ++ ) 
          Rcout << param[offset + i] << " ";
        offset += pW;
      }
      if (zi[2]) {
        Rcout << "), gamma3 = (";
        if (zi[2]) for (int i = 0; i < pW; i ++ ) 
          Rcout << param[offset + i] << " ";
      }
      Rcout << ")" << endl;
    }

    // Updating alpha and eta1 vectors
    // initializing AEnew (alpha and eta1)
    for (int i = 0; i < pZ; i++) {
      AEnew[i] = eta1[i];
    }
    for (int i = 0; i < 3; i++) {
      AEnew[i + pZ] = log(alpha[i]);
    }
    // Rcout << "M-step: alpha, eta1" << endl;
    error = 1.0;
    while ((error >= EPSILON2)) {
      
      
#ifdef DEBUG8
      Rcout << "eta1 before" << endl;
      Rcout << "  eta1 error = "<< error << "/ (e1, a) = ";
      for (int i = 0; i < pZ; i++) {
        Rcout << eta1[i] << " ";
      }
      for (int i = 0; i < 3; i++) {
        Rcout << alpha[i] << " ";
      }
      Rcout << "/ scr = ";
      for (int i = 0; i < pZ + 3; i++) {
        Rcout << score(i, 0) << " ";
      }
      Rcout << endl;
#endif
      
      // update alpha and eta1
      updateAE(ZZ, expt, exptBar, b1, eta1, alpha, 
               AEnew, V_eta1, Zbar, score, n, pZ, error);
      
#ifdef DEBUG8
  Rcout << "  eta1 error = "<< error << "/ (e1, a) = ";
  for (int i = 0; i < pZ; i++) {
    Rcout << eta1[i] << " ";
  }
  for (int i = 0; i < 3; i++) {
    Rcout << alpha[i] << " ";
  }
  Rcout << "/ scr = ";
  for (int i = 0; i < pZ + 3; i++) {
    Rcout << score(i, 0) << " ";
  }
  Rcout << endl;
#endif
  


      // update b1 correspondingly (b1 should be updated within the optimization loop)
      logLin(ZZ, eta1, n, pZ, 1, b1);
    } 

    // Rcout << "M-step: epsilon" << endl;
    error = 1.0;
    while ((error >= EPSILON2)) {
      // update epsilon
// Rcout << "error = "  << error << endl;
Rcout << "V_eps = " << endl;
for (int i = 0; i < pZ; i++) {
  for (int j = 0; j < pZ; j++) {
    Rcout << V_eps(i, j) << " ";
  }
  Rcout << endl;
}
      updateEpsilon(ZZ, expt, b21, epsilon, 
                    epsNew, V_eps, Zbar, score_eps, n, pZ, error);
//#ifdef DEBUG8
Rcout << "   epsilon[0] = " << epsilon[0] << " error = "<< error << endl;    
//#endif
      
      // update eta2 correspondingly
      for (int i = 0; i < pZ; i++) {
        eta2[i] = eta1[i] + epsilon[i];
      }
      // update b21 correspondingly (b21 should be updated within the optimization loop)
      logLin(ZZ, epsilon, n, pZ, 1, b21);
    }

    // After finding the optimal b21, update b2. (b2 need not be updated within the optimization loop)
    for (int i = 0; i < n; i++) {
      b2[i] = b21[i] * b1[i];
    }
    // Rcout << "M-step: gamma" << endl;
    // int counter = 0L;
    error = zi[3] > 0 ? 1.0: 0.0;
    while ((error >= EPSILON3)) {
    // counter++;
    // if (counter > 10000) break;
    // if (0)
    //   {
    //     Rcout << "score_gam = " << endl;
    //     for (int j =0; j < dim_gam; j++) {
    //       Rcout << score_gam[j] << " ";
    //     }
    //     Rcout << endl;
    //     Rcout << "V_gam = " << endl;
    //     for (int i =0; i < dim_gam; i++) {
    //       for (int j =0; j < dim_gam; j++) {
    //         Rcout << V_gam[i + j * dim_gam] << " ";
    //       }
    //       Rcout << endl;
    //     }
    //   }
    
      // update epsilon
      updateGamma(WW, expt, p2, p3, p4, gamma1, gamma2, gamma3, 
                  gammaNew, V_gam, score_gam, n, pW, error, zi);
#ifdef DEBUG8g
  Rcout << "gamma error = ("<< error;
  if (zi[0]) Rcout << ") " << endl << "gamma1 = (";
  if (zi[0]) for (int i =0; i < pW; i++) Rcout << gamma1[i] << ", ";
  if (zi[1]) Rcout << ") " << endl << "gamma2 = ";
  if (zi[1]) for (int i =0; i < pW; i++) Rcout << gamma2[i] << ", ";
  if (zi[2]) Rcout << ") " << endl << "gamma3 = ";
  if (zi[2]) for (int i =0; i < pW; i++) Rcout << gamma3[i] << ", ";
  Rcout << ") " << endl << "score_gam = ";
  if (zi[3] > 0) {
    for (int i = 0; i < dim_gam; i++) {
      Rcout << score_gam(i, 0) << " ";
    }
    Rcout << endl;
  }
  // if (zi[3] > 0) {
  //   Rcout << "  pi2 = ";
  //   for (int i = 0; i < n; i++) {
  //     if (i < 5) Rcout << p2[i] << " (i =" << i <<") " ;
  //   }
  //   Rcout << endl;
  // }
  // Rcout << "(gamma1) = ";
  // for (int i= 0; i < pW; i++) {
  //   Rcout << gamma1[i] << " ";
  // }
  // Rcout << endl;
#endif
      // update pi correspondingly
      if (zi[0]) logLin(WW, gamma1, n, pW, 1, p2);
      if (zi[1]) logLin(WW, gamma2, n, pW, 1, p3);
      if (zi[2]) logLin(WW, gamma3, n, pW, 1, p4);
      for (int i= 0; i < n; i++) { // normalization and p4
        p1[i] = 1.0;
        if (zi[0]) p1[i] += p2[i];
        if (zi[1]) p1[i] += p3[i];
        if (zi[2]) p1[i] += p4[i];
        
        if (zi[0]) p2[i] /= p1[i];
        if (zi[1]) p3[i] /= p1[i];
        if (zi[2]) p4[i] /= p1[i];
        p1[i] = 1.0 / p1[i];
      }
    }
  
    // assigning back each parameter
    int offset = 0;
    for (int i = 0; i < 3; i++) param[offset + i] = alpha[i];
    offset = 3;
    for (int i = 0; i < pZ; i++) param[offset + i] = eta1[i];
    offset += pZ;
    for (int i = 0; i < pZ; i++) param[offset + i] = eta2[i]  ;
    offset += pZ;
    if (zi[0]) {
      for (int i = 0; i < pW; i++) param[offset + i] = gamma1[i];
      offset += pW;
    }
    if (zi[1]) {
      for (int i = 0; i < pW; i++) param[offset + i] = gamma2[i];
      offset += pW;
    }
    if (zi[2])
      for (int i = 0; i < pW; i++) param[offset + i] = gamma3[i];
    
    
    param_diff = 0.0;
    for(int i = 0;i < 3;i++)
    {
      double dif = fabs(0.0 + param[i] - param_old[i]); // for alpha's, original scale
      if( dif > param_diff) param_diff = dif; // maximum difference
    }
    for(int i = 3;i < dim_param;i++)
    {
      double dif = fabs(0.0 + exp(param[i]) - exp(param_old[i])); // for others, exponential scale
      if( dif > param_diff) param_diff = dif; // maximum difference
    }
    
    // updating max lik
    if (exptBar[0] > ll_max) {
      iter_max = iter[0];
      for(int i = 0;i < dim_param; i++) {
        param_max[i] = param_old[i];
      }
      for (int l=0; l < n; l++) {
        for(int i = 0;i < 12;i++) {
          expt_max(i, l) = expt(i, l);
        }
      }
    }
    
    //If iteration lasts more than ITER_ALLOWANCE times without improvement of lik, stops.
    if (iter[0] >= iter_max + ITER_ALLOWANCE) { 
      break;
    }
    trajectory[iter[0]] = expt[0];
    iter[0] += 1;
  }
  
  if (maxiter <= iter[0] &&  param_diff > tol) {nonconv[0] = 1;}  
  
  // cout << "pluging in max_lik iter = " << iter_max << " lik = " << expt_max[0] << " instead of lik = " << expt[0] << endl;
  // replacing back final result with historical maximum
  iter[0] = iter_max;
  for(int i = 0;i < dim_param;i++) {
    param[i] = param_max[i];
  }
  // for (int l = 0; l < n; l++) {
  //   for (int i = 0; i < 12;i++) {
  //     expt[i, l] = expt_max[i, l];
  //   }
  // }

  if (se == 1) { //updating expt and calculate SE when called for.
    dBvZINB_Expt_mat(xvec, yvec, ZZ, WW, n, pZ, pW,
                     alpha, b1, b2, p1, p2, p3, p4, 
                     expt, s_i, s_i_abp, info, se, zi[3], expt_i);
  }

  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
  // for (int i = 0; i < 9; i++) {
  //   param2[i] = param[i];  //returning param to param2
  // }
  
  // List z = List::create(param, xvec, yvec, freq, n, expt, info, se, 
  //                       iter, nonconv, trajectory, zi);
  List z = List::create(param, expt, info, iter, nonconv, trajectory);
  
  return z;
}