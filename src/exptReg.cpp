// [[Rcpp::depends(BH)]]
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <string>
#include "bzinb.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std;
using namespace Rcpp;
using namespace boost::numeric::ublas;

// [[Rcpp::export]]
void logLin(NumericMatrix &mat, NumericVector &vec, int &n, int &p, int sign, 
            NumericVector &result) {
  // mat = n x p matrix
  // vec = p x 1 vector
  // result = exp(sign * mat x vec) = n x 1 vector
  
  for (int i = 0; i < n; i++) {
    result[i] = 0.0;
    for (int j = 0; j < p; j++)
      result[i] += mat(i,j) * vec(j);
    result[i] = exp(sign * result[i]);
  }
}


// [[Rcpp::export]]
void dBvZINB_Expt_mat(IntegerVector &xvec, IntegerVector &yvec,
                      NumericMatrix &ZZ, NumericMatrix &WW,
                      int &n, int &pZ, int &pW,
                      NumericVector &alpha, NumericVector &b1, NumericVector &b2, 
                      NumericVector &p1, NumericVector &p2, 
                      NumericVector &p3, NumericVector &p4,
                      NumericMatrix &expt, NumericVector &s_i, NumericVector &s_i_abp,
                      NumericVector &info, int se, int bnb,
                      NumericVector &expt_i) {
  int x, y, f = 1L;
  double a0 = alpha[0], a1 = alpha[1], a2 = alpha[2];
  double b1i, b2i, p1i, p2i, p3i, p4i;
  // NumericVector expt_i(12);
  
  // initialize result
  for (int i = 0; i < 12 * n; i++) {
    expt[i] = 0.0;
  }
  
  for (int i = 0; i < n; i++) {
    x = xvec[i];
    y = yvec[i];
    b1i = b1[i];
    b2i = b2[i];
    p1i = p1[i];
    p2i = p2[i];
    p3i = p3[i];
    p4i = p4[i];
    // expt_i = &expt + i * 12;
    expt_i = expt( _, i);
    // sumFreq += freq[i];

// Rcout << "before " << endl;
// for (int j=0; j<12; j++) {
//   Rcout << "j = " << j << ", expt = " << expt[j + i * 12] << " " << endl;  
// }
    
    // add expectations with weights (freq)
    dBvZINB_Expt(x, y, f, a0, a1, a2,
                 b1i, b2i, p1i, p2i, p3i, p4i,
                 // b1[i], b2[i], p1[i], p2[i], p3[i], p4[i],
                 // &b1 + i, &b2 + i, &p1 + i, &p2 + i, &p3 + i, &p4 + i,
                 expt_i, s_i_abp, info, se, bnb, 1);
                 //*expt_i, s_i, info, se, bnb);
    // expt( _, i) = expt_i;
    for (int j = 0; j < 12; j++) expt(j, i) = expt_i[j];
    
    
    if (se) {
      int offset;
      for (int j = 0; j < 3; j++) s_i[j] = s_i_abp[j];
      offset = 3;
      for (int j = 0; j < pZ; j++) s_i[j + offset] = s_i_abp[3] * b1i * ZZ[i + n * j];
      offset += pZ;
      //for (int j = 0; j < pZ; j++) s_i[j + offset] = s_i_abp[4] * b2i/b1i * ZZ[i + n * j]; // for epsilon
      for (int j = 0; j < pZ; j++) s_i[j + offset] = s_i_abp[4] * b2i * ZZ[i + n * j]; // for eta2
      offset += pZ;
      int offset2 = offset + pW;
      int offset3 = offset2 + pW;
      for (int j = 0; j < pW; j++) {
          s_i[j + offset]   = s_i_abp[5] * p1i * (1-p1i) * WW[i + n * j];
          s_i[j + offset]  += s_i_abp[6] * p1i * ( -p2i) * WW[i + n * j];
          s_i[j + offset]  += s_i_abp[7] * p1i * ( -p3i) * WW[i + n * j];
          s_i[j + offset2]  = s_i_abp[5] * p2i * ( -p1i) * WW[i + n * j];
          s_i[j + offset2] += s_i_abp[6] * p2i * (1-p2i) * WW[i + n * j];
          s_i[j + offset2] += s_i_abp[7] * p2i * ( -p3i) * WW[i + n * j];
          s_i[j + offset3]  = s_i_abp[5] * p3i * ( -p1i) * WW[i + n * j];
          s_i[j + offset3] += s_i_abp[6] * p3i * ( -p2i) * WW[i + n * j];
          s_i[j + offset3] += s_i_abp[7] * p3i * (1-p3i) * WW[i + n * j];
      }
      
      offset = offset3 + pW;
      for (int j = 0; j < offset; j++) {
        for (int k = 0; k < offset; k++) {
          info[j + offset*k] += s_i[j] * s_i[k];
        }
      }
      
//info[0,0] = s_i[0]
    }
    // if (i==0) Rcout << "expt" << endl;
    // for (int j=0; j<12; j++) {
    //   Rcout << expt_i[j] << " ";
    // }
    // Rcout << endl;
    
// Rcout << "after" << endl;
// for (int j=0; j<12; j++) {
//   Rcout << "j = " << j << ", expt = " << expt[j + i * 12] << " " << endl;  
// }
// Rcout << "after (expt_i)" << endl;
// for (int j=0; j<12; j++) {
//   Rcout << "j = " << j << ", expt = " << expt_i[j] << " " << endl;  
// }

  }
}
