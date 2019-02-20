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

#define EPSILON 1e-7
#define EPSILON2 1e-7

// 2. optimization
double inv_digamma(double x, double y) 
{ 
  double h = (boost::math::digamma(x) - y) / boost::math::trigamma(x);
  while (fabs(h) >= EPSILON2)
  {
    h = (boost::math::digamma(x) - y) / boost::math::trigamma(x);
    while (h > x) {
      // cout <<"inv_digamm hit zero. ";
      h /= 2.0;
    }
    x = x - h;
  }
  //cout << endl << "inv_digamm = " << x << endl;
  return(x);
  //cout << "The value of x is : " << x <<;
}


// optimization


void inv_digamma_vec(NumericVector& lb, NumericVector &expt, NumericVector &a, NumericVector &idgam) 
{ 
  //double idgam[3];
  //double *idgam = new double[3];
  idgam[0] = inv_digamma(a[0], expt[4] - lb[0]);
  idgam[1] = inv_digamma(a[1], expt[5] - lb[0]);
  idgam[2] = inv_digamma(a[2], expt[6] - lb[0]);
  //return(idgam);
}

// double objective(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam) 
// { 
//   inv_digamma_vec(lb, expt, a, idgam);
//   // cout << "idgam123: "<< idgam[0] << " "<< idgam[1] << " "<< idgam[2] << endl;
//   //  inv_digamma_vec(b1, expt, a, idgam);
//   // double result = (double);
//   return (log(idgam[0] + idgam[1] + idgam[2]) + lb - log(expt[1] + expt[2] + expt[3]));
// }
// 
// // Derivative
// double derivFunc(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam)
// {
//   inv_digamma_vec(lb, expt, a, idgam);
//   double result = 0.0;
//   for (int i = 0; i < 3; i++) {
//     result += 1/ boost::math::trigamma(idgam[i]);
//   }
//   return ((result)/(idgam[0] + idgam[1] + idgam[2]) + 1.0);
// }

// Derivative
double hfunc(NumericVector& lb, NumericVector &expt, NumericVector &a, NumericVector &idgam)
{
  inv_digamma_vec(lb, expt, a, idgam);
  double result = 0.0;
  for (int i = 0; i < 3; i++) {
    // cout << "idgam[0:2]" << idgam[0] << " " << idgam[1] << " " << idgam[2] << endl;
    result += (1/ boost::math::trigamma(idgam[i]));
  }
  result = - result / (idgam[0] + idgam[1] + idgam[2]) + 1.0;
  result = (log(idgam[0] + idgam[1] + idgam[2]) + lb[0] - log(expt[1] + expt[2] + expt[3])) / result;
  return (result);
}

// [[Rcpp::export]]
void opt_lb(NumericVector& lb, NumericVector &expt, NumericVector &a,  NumericVector &idgam)
{
  //double* lb = log(b1);
  // double h =  objective(lb, expt, a, idgam) / derivFunc(lb, expt, a, idgam);
  double h = hfunc(lb, expt, a, idgam);
  // cout << "h = " << h << " ";
  while (fabs(h) >= EPSILON)
  {
    // h = objective(lb, expt, a, idgam)/derivFunc(lb, expt, a, idgam);
    // while (h > b1) {
    //   // cout << "opt_b1 hit zero. ";
    //   h /= 2.0;
    // }
    // cout << "h = " << h << " ";
    lb[0] -= h;
    h = hfunc(lb, expt, a, idgam);
  }
  //b1 = exp(lb);
  //cout << "lb =" << lb  << ", b1 = " << exp(lb) << endl;
  // double* idgam = inv_digamma_vec(lb, expt, a);
  a[0] = idgam[0];
  a[1] = idgam[1];
  a[2] = idgam[2];
  // cout << "a123: "<< a[0] << " "<< a[1] << " "<< a[2] << " " << exp(lb) << endl;
}