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
                      NumericMatrix &expt, NumericVector &s_i, 
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
                 expt_i, s_i, info, se, bnb);
                 //*expt_i, s_i, info, se, bnb);
    // expt( _, i) = expt_i;
    for (int j = 0; j < 12; j++) expt(j, i) = expt_i[j];
    
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

// if (i>0) break;
  }
}
