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
// #define DEBUG4
// [[Rcpp::depends(BH)]]

// 3. EM
// [[Rcpp::export]]
void em(NumericVector& param, const IntegerVector &xvec, const IntegerVector &yvec, 
        const IntegerVector &freq, const int &n, NumericVector &expt, NumericVector &info,
        const int &se, IntegerVector &iter, int &maxiter, double &tol, int showFlag)
{
  double param_diff = 1.0;
  NumericVector param_old(9);
  NumericVector idgam(3);
  NumericVector lb(1);
  NumericVector s_i(8, 0.0);
  iter[0] = 0;
  
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

    iter[0] += 1;
    //cout << "maxiter = " << maxiter << " iter1 = " << iter[0] << endl;  
    for(int i = 0;i < 9;i++)
    {
      param_old[i] = param[i];
    }
    
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, 0);
    double delta = expt[11]*1.0 / (expt[1] + expt[3]);
    
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
      //cout <<  i <<": "<< param[i] <<  ", dif = " << dif <<endl;
      //double dif = fabs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff)
      {
        param_diff = dif;
      }
    }
  }
  
  
  if (se == 1) { //updating expt and calculate SE when called for.
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, se);
  }
  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
}