#include <Rcpp.h>
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
// #define DEBUG4
// [[Rcpp::depends(BH)]]

// 3. EM
// [[Rcpp::export]]
void em(NumericVector& param2, const IntegerVector &xvec, const IntegerVector &yvec, 
        const IntegerVector &freq, const int &n, NumericVector &expt, NumericVector &info,
        const int &se, IntegerVector &iter, int &maxiter, double &tol, int showFlag,
        NumericVector& trajectory)
{
  long double param[9];
  double param_diff = 1.0;
  long double param_old[9];
  long double idgam[3];
  long double lb[1];
  NumericVector s_i(8, 0.0);
  iter[0] = 0;
  
  // storage for max_likelihood.
  double param_max[9];
  int iter_max = 0;
  double expt_max[12];
  
  for (int i = 0; i < 9; i++) {
    param[i] = (long double) param2[i];  //initializing param2 with param
  }
  
  //cout << "maxiter = " << maxiter << " iter = " << iter[0] << endl;  
  while(maxiter > iter[0] && param_diff > tol)
  {
    // cout << "param_diff = " << param_diff << " ";
    if (showFlag > 0 & iter[0] > showFlag) {
      Rcout << "iter = " << iter[0]  << ", likelihood = " << expt[0] << endl <<
      "  (a0-2, b1-2, p1-4) = (" << param[0] << ", " << param[1] << ", " << param[2] << 
        ", " << param[3] << ", " << param[4] << ", " << param[5] << ", " << param[6] << 
          ", " << param[7] << ", " << param[8] << ")" << endl;
    }
    if (iter[0] < 3) {
      expt_max[0] = expt[0];
    }
    //cout << "maxiter = " << maxiter << " iter1 = " << iter[0] << endl;  
    for(int i = 0;i < 9;i++)
    {
      param_old[i] = param[i];
      // cout << "param [" << i << "] = " << param[i] << " ";
    }
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, 0);
    
    long double delta = expt[11]*1.0 / (expt[1] + expt[3]);
    
    for(int i = 0;i < 4; i++)
    {
      param[5 + i] = expt[7 + i];
    }
    lb[0] = log(param[3]);
#ifdef DEBUG4
    for (int i = 0; i < 12; i++) 
    {
      if (i == 0) {Rcout << "expt: ";}
      Rcout << expt[i] << " ";
      if (i == 11) {Rcout << "lb = " << lb[0] << endl;}
    }
    for (int i = 0; i < 8; i++) 
    {
      if (i == 0) {Rcout << "param: ";}
      Rcout << param[i] << " ";
      if (i == 7) {Rcout <<  endl;}
    }
    for (int i = 0; i < 3; i++) 
    {
      if (i == 0) {Rcout << "idgam: ";}
      Rcout << idgam[i] << " ";
      if (i == 3) {Rcout <<  endl;}
    }
    Rcout << "before opt_lb! (of iter" << iter << ") ";
#endif
    //if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    // Finding optimized a0, a1, a2, b1
    opt_lb(lb, expt, param, idgam);
#ifdef DEBUG4
  Rcout << ", after opt_lb! (of iter" << iter << ") "<< endl;
#endif
    //if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    param[3] = exp(lb[0]);
    param[4] = param[3] * delta;
    
    param_diff = 0.0;
    
    for(int i = 0;i < 9;i++)
    {
      double dif = fabs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff) param_diff = dif;
    }
    
    // updating max lik
    if (expt[0] > expt_max[0]) {
      iter_max = iter[0];
      for(int i = 0;i < 9;i++) {
        param_max[i] = param_old[i];
      }
      for(int i = 0;i < 12;i++) {
        expt_max[i] = expt[i];
      }
    }
    
    //If iteration lasts more than ITER_ALLOWANCE times without improvement of lik, stops.
    if (iter[0] >= iter_max + ITER_ALLOWANCE) { 
      break;
    }
    trajectory[iter[0]] = expt[0];
    iter[0] += 1;
  }
  
  // cout << "pluging in max_lik iter = " << iter_max << " lik = " << expt_max[0] << " instead of lik = " << expt[0] << endl;
  // replacing back final result with historical maximum
  iter[0] = iter_max;
  for(int i = 0;i < 9;i++) {
    param[i] = param_max[i];
  }
  for(int i = 0;i < 12;i++) {
    expt[i] = expt_max[i];
  }

  if (se == 1) { //updating expt and calculate SE when called for.
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, se);
  }

  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
  for (int i = 0; i < 9; i++) {
    param2[i] = param[i];  //returning param to param2
  }
}