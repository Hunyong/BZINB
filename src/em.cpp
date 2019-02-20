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
// [[Rcpp::depends(BH)]]

// 3. EM
// [[Rcpp::export]]
void em(NumericVector& param, const IntegerVector &xvec, const IntegerVector &yvec, 
        const IntegerVector &freq, const int &n, NumericVector &expt, 
        IntegerVector &iter, int &maxiter, double &tol, int showFlag)
{
  double param_diff = 1.0;
  NumericVector param_old(9);
  NumericVector idgam(3);
  NumericVector lb(1);
  iter[0] = 0;
  
  //cout << "maxiter = " << maxiter << " iter = " << iter[0] << endl;  
  while(maxiter > iter[0] && param_diff > tol)
  {
    // cout << "param_diff = " << param_diff << " ";
    if (showFlag == 1 && iter[0] >2059) {
      cout << "iter = " << iter[0]  << ", a0 = " << param[0] << ", a1 = " << param[1]
           << ", a2 = " << param[2] << ", b1 = " << param[3] << ", b2 = " << param[4] 
           << ", p1 = " << param[5] << ", p2 = " << param[6] << ", p3 = " <<  param[7]
           << ", p4 = " << param[8] << ", likelihood = " << expt[0] << endl;
    }
    iter[0] += 1;
    //cout << "maxiter = " << maxiter << " iter1 = " << iter[0] << endl;  
    for(int i = 0;i < 9;i++)
    {
      param_old[i] = param[i];
    }
    
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt);
    double delta = expt[11]*1.0 / (expt[1] + expt[3]);
    
    for(int i = 0;i < 4; i++)
    {
      param[5 + i] = expt[7 + i];
    }
    lb[0] = log(param[3]);
    if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    // Finding optimized a0, a1, a2, b1
    opt_lb(lb, expt, param, idgam);
    if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    param[3] = exp(lb[0]);
    param[4] = param[3] * delta;
    
    param_diff = 0.0;
    for(int i = 0;i < 9;i++)
    {
      double dif = fabs(0.0 + param[i] - param_old[i]);
      cout <<  i <<": "<< param[i] <<  ", dif = " << dif <<endl;
      //double dif = fabs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff)
      {
        param_diff = dif;
      }
    }
  }
  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
}