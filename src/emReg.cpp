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
#define EPSILON3 1e-5  //for M-steps (gamma)
// #define DEBUG4
// #define DEBUG5
// #define DEBUG7
// #define DEBUG8

// 3. EM
// [[Rcpp::export]]
List emReg(NumericVector& param2, IntegerVector &xvec, IntegerVector &yvec, 
           NumericMatrix& ZZ, NumericMatrix& WW,
           int &pZ, int &pW,
           int &n, int &se, int &maxiter, double &tol, int showFlag,
           int bnb)
{
  NumericVector param = clone(param2);
  int dim_param = 3 + 2 * pZ + 3 * pW;
  NumericVector alpha(3);
  NumericVector eta1(pZ);
  NumericVector eta2(pZ);
  NumericVector gamma1(pW);
  NumericVector gamma2(pW);
  NumericVector gamma3(pW);
  double error;

  double param_diff = 1.0;
  NumericVector param_old(dim_param);
  //double idgam[3];
  //double lb[1];  // replaced by eta1
  NumericMatrix expt(12, n);
  NumericVector exptBar(12);
  NumericVector b1(n);
  NumericVector b2(n);
  NumericVector p1(n);
  NumericVector p2(n);
  NumericVector p3(n);
  NumericVector p4(n);
  NumericVector b21(n); // b2 / b1
  NumericVector Zbar(pZ);
  NumericVector epsilon(pZ);
  NumericVector gamma(pW * 3);
  NumericVector gammaNew(pW * 3);
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
  arma::mat V_gam(pW * 3, pW * 3);
  arma::mat score_gam(pW * 3, 1);
  
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
  for (int i = 0; i < pW; i++) gamma1[i] = (double) param[offset + i];
  offset += pW;
  for (int i = 0; i < pW; i++) gamma2[i] = (double) param[offset + i];
  offset += pW;
  for (int i = 0; i < pW; i++) gamma3[i] = (double) param[offset + i];
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
    for(int i = 0;i < dim_param;i++)
    {
      param_old[i] = param[i];
      // cout << "param [" << i << "] = " << param[i] << " ";
    }

    if (iter[0] <= 1L) {
      // reparametrization
      logLin(ZZ, eta1, n, pZ, 1, b1);
      logLin(ZZ, eta2, n, pZ, 1, b2);
      for (int i = 0; i < n; i++) b21[i] = b2[i] / b1[i];
      logLin(WW, gamma1, n, pW, 1, p1);
      logLin(WW, gamma2, n, pW, 1, p2);
      logLin(WW, gamma3, n, pW, 1, p3);
      for (int i= 0; i < n; i++) { // normalization and p4
        p4[i] = p1[i] + p2[i] + p3[i] + 1.0;
        p1[i] = p1[i] / p4[i];
        p2[i] = p2[i] / p4[i];
        p3[i] = p3[i] / p4[i];
        p4[i] = 1.0 / p4[i];
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

    dBvZINB_Expt_mat(xvec, yvec, ZZ, WW, n, pZ, pW,
                     alpha, b1, b2, p1, p2, p3, p4, 
                     expt, s_i, s_i_abp, info, 0, bnb, expt_i);

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
      Rcout << "), gamma1 = (";
      for (int i = 0; i < pW; i ++ ) 
        Rcout << param[3 + 2 * pZ + i] << " ";
      Rcout << "), gamma2 = (";
      for (int i = 0; i < pW; i ++ ) 
        Rcout << param[3 + 2 * pZ + pW + i] << " ";
      Rcout << "), gamma3 = (";
      for (int i = 0; i < pW; i ++ ) 
        Rcout << param[3 + 2 * pZ + 2 * pW + i] << " ";
      Rcout << ")";
      // 
      // if (!bnb) {
      //   Rcout << ", p1 " << param[5] << ", p2 " << param[6] << ", p3 " << param[7] <<
      //     ", p4 " << param[8];
      // }
      Rcout << endl;
    }

    // #ifdef DEBUG4
    //     for (int i = 0; i < 8; i++) 
    //     {
    //       if (i == 0) {Rcout << "param: ";}
    //       Rcout << param[i] << " ";
    //       if (i == 7) {Rcout <<  endl;}
    //     }
    //     for (int i = 0; i < 12; i++) 
    //     {
    //       if (i == 0) {Rcout << "expt: ";}
    //       Rcout << expt[i] << " ";
    //     }
    //     Rcout << endl;
    // #endif
    
    /// HERE...opt_...!!!!!
    // Updating alpha and eta1 vectors
    // initializing AEnew (alpha and eta1)
    for (int i = 0; i < pZ; i++) {
      AEnew[i] = eta1[i];
    }
    for (int i = 0; i < 3; i++) {
      AEnew[i + pZ] = log(alpha[i]);
    }

    error = 1.0;
    while ((error >= EPSILON2)) {
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


    error = 1.0;
    while ((error >= EPSILON2)) {
      // update epsilon
// Rcout << "error = "  << error << endl;
      updateEpsilon(ZZ, expt, b21, epsilon, 
                    epsNew, V_eps, Zbar, score_eps, n, pZ, error);
#ifdef DEBUG8
  Rcout << "   epsilon[0] = " << epsilon[0] << " error = "<< error << endl;    
#endif
      
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
    
// int counter = 0L;
    error = 1.0;
    while ((error >= EPSILON3)) {
// counter++;
// if (counter > 500) break;
// if (counter >= 490) {
//   Rcout << "score_gam = " << endl;
//   for (int j =0; j < 3*pW; j++) {
//     Rcout << score_gam[j] << " ";
//   }  
//   Rcout << endl;
//   Rcout << "V_gam = " << endl;
//   for (int i =0; i < 3*pW; i++) {
//     for (int j =0; j < 3*pW; j++) {
//       Rcout << V_gam[i + j * 3*pW] << " ";
//     }  
//     Rcout << endl;
//   }
// }
      // update epsilon
      updateGamma(WW, expt, p1, p2, p3, gamma1, gamma2, gamma3, 
                  gammaNew, V_gam, score_gam, n, pW, error);
#ifdef DEBUG8
  Rcout << "gamma error = "<< error <<" ";
  for (int i =0; i < pW; i++) Rcout << "("<< gamma1[i] << ", "<< gamma2[i] << ", "<< gamma3[i] << ") ";
  Rcout << "/ scr = ";
  for (int i = 0; i < pW* 3; i++) {
    Rcout << score(i, 0) << " ";
  }
  Rcout << endl;
#endif
      // update pi correspondingly
      logLin(WW, gamma1, n, pW, 1, p1);
      logLin(WW, gamma2, n, pW, 1, p2);
      logLin(WW, gamma3, n, pW, 1, p3);
      for (int i= 0; i < n; i++) { // normalization and p4
        p4[i] = p1[i] + p2[i] + p3[i] + 1.0;
        p1[i] = p1[i] / p4[i];
        p2[i] = p2[i] / p4[i];
        p3[i] = p3[i] / p4[i];
        p4[i] = 1.0 / p4[i];
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
    for (int i = 0; i < pW; i++) param[offset + i] = gamma1[i];
    offset += pW;
    for (int i = 0; i < pW; i++) param[offset + i] = gamma2[i];
    offset += pW;
    for (int i = 0; i < pW; i++) param[offset + i] = gamma3[i];
    
    param_diff = 0.0;
    for(int i = 0;i < dim_param;i++)
    {
      double dif = fabs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff) param_diff = dif;
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
                     expt, s_i, s_i_abp, info, se, bnb, expt_i);
  }

  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
  // for (int i = 0; i < 9; i++) {
  //   param2[i] = param[i];  //returning param to param2
  // }
  
  // List z = List::create(param, xvec, yvec, freq, n, expt, info, se, 
  //                       iter, nonconv, trajectory, bnb);
  List z = List::create(param, expt, info, iter, nonconv, trajectory);
  
  return z;
}