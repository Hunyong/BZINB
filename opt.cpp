// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#define EPSILON 1e-6

// using namespace std;
// using namespace Rcpp;


double func(double x, double y) 
{ 
  return (boost::math::digamma(x) - y);
}

double inv_digamma(double x, double y) 
{ 
  double h = func(x, y) / boost::math::trigamma(x);
  while (abs(h) >= EPSILON)
  {
    h = func(x, y) / boost::math::trigamma(x);
    x = x - h;
  }
  return(x);
  //cout << "The value of x is : " << x;
}


// optimization


double* inv_digamma_vec(double b1, NumericVector &expt, NumericVector &a) 
{ 
  static double idgam[3];
  //double *idgam = new double[3];
  idgam[0] = inv_digamma(a[0], expt[3] + log(b1));
  idgam[1] = inv_digamma(a[1], expt[4] + log(b1));
  idgam[2] = inv_digamma(a[2], expt[5] + log(b1));
  return(idgam);
}

double objective(double b1, NumericVector &expt, NumericVector &a) 
{ 
  double* idgam = inv_digamma_vec(b1, expt, a);
  // cout << "idgam123: "<< idgam[0] << " "<< idgam[1] << " "<< idgam[2] << endl;
  //  inv_digamma_vec(b1, expt, a, idgam);
  // double result = (double);
  return ((idgam[0] + idgam[1] + idgam[2]) * b1 - expt[0] - expt[1] - expt[2]);
}

// Derivative
double derivFunc(double b1, NumericVector &expt, NumericVector &a)
{
  double* result = inv_digamma_vec(b1, expt, a);
  result[0] += 1/ boost::math::trigamma(result[0]);
  result[1] += 1/ boost::math::trigamma(result[1]);
  result[2] += 1/ boost::math::trigamma(result[2]);
  return (result[0] + result[1] + result[2]);
}

// [[Rcpp::export]]
void opt_b1(double& b1, NumericVector &expt, NumericVector &a)
{
  double h = objective(b1, expt, a) / derivFunc(b1, expt, a);
  while (abs(h) >= EPSILON)
  {
    h = objective(b1, expt, a)/derivFunc(b1, expt, a);
    b1 = b1 - h;
  }
  double* idgam = inv_digamma_vec(b1, expt, a);
  a[0] = idgam[0];
  a[1] = idgam[1];
  a[2] = idgam[2];
  //cout << "a123: "<< a[0] << " "<< a[1] << " "<< a[2] << " " << b1 << endl;
}