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

// void e_Z_eta(NumericMatrix &ZZ, NumericVector &eta1, NumericVector &result, int sign, int n, int p) {
//   // ZZ = n x p matrix
//   // eta1 = p x 1 vector
//   // result = exp(sign * ZZ eta1)
//   
//   for (int i = 0; i < n; i++) {
//     
//   }
//   
// }

// update alpha and eta1
void updateAE(NumericMatrix &ZZ, NumericMatrix &expt,  NumericVector &exptBar, NumericVector &b1, 
              NumericVector &eta1, NumericVector &alpha, NumericVector &AEnew, 
              arma::mat &V_eta1, NumericVector &Zbar, arma::mat &score, 
              int &n, int &p, double &error) {
  // ZZ = n x p matrix
  // expt = 12 x n matrix
  // eta1 = pZ x 1 vector
  // alpha = 3 x 1 vector
  // score = (p + 3) x 1 score vector
  //       = Z' (-D (a0 + a1 + a2) + D_{R0 + R1 + R2} D_{exp(-Z eta1)}) 1n   | (D_{R0 - Z eta1} 1n - digamma(a0) |  (D_{R1 - Z eta1} 1n - digamma(a1) |  (D_{R2 - Z eta1} 1n - digamma(a2) 
  // V_eta1 = (p + 3) x (p + 3) Hessian
  //        = (Z' D_{R0 + R1 + R2} D_{exp(-Z eta1)} Z | n Z_bar' | n Z_bar' | n Z_bar')
  //           n Z_bar                                | n trigamma(a0) | 0 | 0
  //           n Z_bar                                | 0 | n trigamma(a1) | 0
  //           n Z_bar                                | 0 | 0 | n trigamma(a2)
  // 
  double R123b1, a012, ZbarEta1 = 0.0, da[3];
  a012 = alpha[0] + alpha[1] + alpha[2];
  for (int i = 0; i < 3; i++) da[i] = boost::math::digamma(alpha[i]);
    
  for (int j = 0; j < p; j++) ZbarEta1 += Zbar[j] * eta1[j]; // = Zbar' eta1
// Rcout << "a012 = " << alpha[0] << " " << alpha[1] << " " << alpha[2] << " " << a012 << endl;  
  
  // initialize V_eta1 and score
  V_eta1.fill(0.0);
  score.fill(0.0);

// Rcout << "p = " << p << endl;  
  // V_eta1 and score: 1 ~ p
  for (int l = 0; l < n; l++) {
    R123b1 = (expt(1, l) + expt(2, l) + expt(3, l))/b1[l];
    for (int i = 0; i < p; i++) {
      score(i,0) += ZZ[l + i * n] * (R123b1 - a012);
      for (int j = 0; j < p; j++) {
        V_eta1(i, j) += ZZ[l + i * n] * R123b1 * ZZ[l + j * n];
      }
    }
  }
// Rcout << "score[0:7] = " << score[0] << " " << score[1] << " " << score[2] << " " << score[3] << " " << endl;
  // V_eta1 and score: p + 1 ~ p + 3
  for (int i = 0; i < 3; i++) {
    score(p + i, 0) = (exptBar[i + 4] - ZbarEta1 - da[i]) * alpha[i];
    for (int j = 0; j < p; j++) {
      V_eta1(p + i, j) = Zbar[j] * alpha[i];
      V_eta1(j, p + i) = Zbar[j] * alpha[i];
    }
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        V_eta1(p + i, p + i) = boost::math::trigamma(alpha[i]) * alpha[i] * alpha[i];
        V_eta1(p + i, p + i) -= score(p + i, 0);
      } else {
        V_eta1(p + i, p + j) = 0.0;
        V_eta1(p + j, p + i) = 0.0;
      }
    }
  }
// Rcout << "V_eta1 = "<< endl;
//   for (int i = 0; i < p + 3; i++) {
//     for (int j = 0; j < p + 3; j++) {
//       Rcout << V_eta1(i, j) << " ";
//     }
//     Rcout << endl;
//   }
  // inverse of V_eta
  V_eta1 = inv(V_eta1);
// Rcout << "inv(V_eta1)[0:2] = "<< V_eta1[0] << " "<< V_eta1[1] << " "<<  V_eta1[2] << endl;

  // updating
  for (int i = 0; i < p + 3; i++) {
    for (int j = 0; j < p + 3; j++) {
      AEnew[i] += V_eta1(i, j) * score(j, 0);
    }
  }
// Rcout << "AEnew = ";
// Rcout << AEnew[0] << " ";

  // getting errors and updating
  error = 0.0;
  for (int i = 0; i < p; i++) {
    error += fabs(AEnew[i] - eta1[i]);
    eta1[i] = AEnew[i];
  }
  for (int i = 0; i < 3; i++) {
    error += fabs(exp(AEnew[i + p]) - alpha[i]);
    alpha[i] = exp(AEnew[i + p]);
  }
}


void updateEpsilon(NumericMatrix &ZZ, NumericMatrix &expt, NumericVector &b21, 
                   NumericVector &epsilon, NumericVector &epsNew, 
                   arma::mat &V_eps, NumericVector &Zbar, arma::mat &score, 
                   int &n, int &p, double &error) {
  // ZZ = n x p matrix
  // expt = 12 x n matrix
  // epsilon = pZ x 1 vector
  // score = (p) x 1 score vector
  //       = Z' (D (a0 + a1 + a2) - D_{R0 + R1 + R2} D_{exp(-Z eta1)}) 1n   | (D_{R0 - Z eta1} 1n - a0 |  (D_{R1 - Z eta1} 1n - a1 |  (D_{R2 - Z eta1} 1n - a2 
  // V_eps = (p) x (p) Hessian inverse
  //        = [  Z' D_{R0 + R2} D_{exp(Z epsilon)} Z  ]^-1
  // 
  double R02eps;
  
  // initialize V_eps and score
  V_eps.fill(0.0);
  score.fill(0.0);
  
  // for (int i = 0; i < p * p; i++) V_eps[i] = 0.0;
  // for (int i = 0; i < p; i++) score[i] = 0.0;
  
  // V_eta1 and score: 1 ~ p
  for (int i = 0; i < p; i++) {
    for (int l = 0; l < n; l++) {
      R02eps = (expt(1, l) + expt(3, l)) * b21[l];
      score(i,0) += ZZ[l + i * n] * (-R02eps + expt(11, l));
      for (int j = 0; j < p; j++) {
        V_eps(i, j) += ZZ[l + i * n] * R02eps * ZZ[l + j * n];
      }
    }
  }
  // inverse of V_eps
  V_eps = inv(V_eps);

  // updating
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      epsNew[i] += V_eps(i, j) * score(j, 0);
    }
  }
  // getting errors and updating
  error = 0.0;
  for (int i = 0; i < p; i++) {
    error += fabs(epsNew[i] - epsilon[i]);
    epsilon[i] = epsNew[i];
  }
}


void updateGamma(NumericMatrix &WW, NumericMatrix &expt,
                 NumericVector &p1, NumericVector &p2, NumericVector &p3, 
                 NumericVector &gamma1, NumericVector &gamma2, NumericVector &gamma3, 
                 NumericVector &gammaNew, 
                 arma::mat &V, arma::mat &score, 
                 int &n, int &p, double &error) {
  // WW = n x p matrix
  // expt = 12 x n matrix
  // gamma = 3p x 1 vector
  // score = (3p) x 1 score vector
  //       = [W' (D_{E1} - D_pi_1) 1]
  //         [W' (D_{E2} - D_pi_2) 1]
  //         [W' (D_{E3} - D_pi_3) 1]
  // V_gam = (3p) x (3p) Hessian inverse
  //        =[ W' [D_pi_1 - D_pi_1^2] W   |  W' [    - D_{pi_1 pi_2}] W  |  W' [   - D_{pi_1 pi_3}] W
  //           W' [  - D_{pi_2 pi_1}] W   |  W' [D_pi_2 - D_{pi_2^2}] W  |  W' [   - D_{pi_2 pi_3}] W
  //           W' [  - D_{pi_3 pi_1}] W   |  W' [    - D_{pi_3 pi_2}] W  |  W' [D_pi_3- D_{pi_3^2}] W ]^-1
  double tmp;
  
  // initialize V and score
  V.fill(0.0);
  score.fill(0.0);
  
  // V_eta1 and score: 1 ~ p
  for (int i = 0; i < p; i++) {
    for (int l = 0; l < n; l++) {
      score(i,       0) += WW[l + i * n] * (expt(7, l) - p1[l]);
      score(i + p,   0) += WW[l + i * n] * (expt(8, l) - p2[l]);
      score(i + 2*p, 0) += WW[l + i * n] * (expt(9, l) - p3[l]);
    }
  }
  for (int i = 0; i < p; i++) {
    for (int l = 0; l < n; l++) {
      for (int j = 0; j < p; j++) {
        V(i, j)             += WW[l + i * n] * (p1[l] - p1[l]*p1[l]) * WW[l + j * n];
        V(i + p, j + p)     += WW[l + i * n] * (p2[l] - p2[l]*p2[l]) * WW[l + j * n];
        V(i + 2*p, j + 2*p) += WW[l + i * n] * (p3[l] - p3[l]*p3[l]) * WW[l + j * n];
        tmp = - WW[l + i * n] * (p1[l] * p2[l]) * WW[l + j * n];
        V(i, j + p)         += tmp;
        V(i + p, j)         += tmp;
        tmp = - WW[l + i * n] * (p1[l] * p3[l]) * WW[l + j * n];
        V(i, j + 2*p)         += tmp;
        V(i + 2*p, j)         += tmp;
        tmp = - WW[l + i * n] * (p2[l] * p3[l]) * WW[l + j * n];
        V(i + p, j + 2*p)         += tmp;
        V(i + 2*p, j + p)         += tmp;
      }
    }
  }
  // inverse of V
  V = inv(V);
  
  // updating
  for (int i = 0; i < 3*p; i++) {
    for (int j = 0; j < 3*p; j++) {
      gammaNew[i] += V(i, j) * score(j, 0);
    }
  }
  // getting errors and updating
  error = 0.0;
  for (int i = 0; i < p; i++) {
    error += fabs(gammaNew[i] - gamma1[i]);
    error += fabs(gammaNew[i + p] - gamma2[i]);
    error += fabs(gammaNew[i + 2*p] - gamma3[i]);
    gamma1[i] = gammaNew[i];
    gamma2[i] = gammaNew[i + p];
    gamma3[i] = gammaNew[i + 2*p];
  }
}