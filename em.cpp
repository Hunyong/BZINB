// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include "exp.cpp"
#include "opt.cpp"

using namespace std;
using namespace Rcpp;

// double dabs(double x) {
//   if (x < 0.0) {
//     x = -x;
//   }
//   return (x);
// }

// [[Rcpp::export]]
void em(NumericVector& param, const IntegerVector &xvec, const IntegerVector &yvec, 
        const IntegerVector &freq, int &n, NumericVector &expt, 
        int &iter,
        int &maxiter, double &tol)
{
  //int iter = 0;
  double param_diff = 10;
  NumericVector param_old[9];
  
  while(maxiter > iter && param_diff > tol)
  {
    iter = iter + 1;
    for(int i = 0;i < 9;i++)
    {
      param_old[i] = param[i];
    }
    
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt);//need to call it first 
    double delta = expt[11]*1.0 / (expt[1] + expt[3]);
    for(int i = 0;i <= 3; i++)
    {
      param[6 + i] = expt[8 + i];
    }
    
    // Finding optimized a0, a1, a2, b1
    opt_b1(param[3], expt, param);
    param[4] = param[3] * delta;
    
    double param_diff;
    for(int i = 0;i < 9;i++)
    {
      double dif = 0;
      cout << param[i] << endl;
      cout << &param << endl;
      //double dif = abs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff)
      {
        param_diff = dif;
      }
    }
  }
}