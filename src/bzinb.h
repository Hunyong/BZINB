#ifndef BZINB_H
#define BZINB_H

#include <Rcpp.h>
#include <math.h>
#include <string>
#include <iostream>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
// 
// double EPSILON1;
// double EPSILON2;

void l1(int& x, int& y, double& a0, double& a1, double& a2, int &k, int& m, 
        long double& result, double adjj = 0);
void l1_c (double& t1, double& t2, int& k, int& m, long double& result, double adjj);
void l1_AC (double& t1, double& t2, int& x, int& y, double& a0, double& a1, double& a2, 
            int& k, int& m, long double& result,  double adjj = 0);
void l2_A (int& x, double& a0, double& a1, double& a2, int& k, long double& result, double adjj);
void l3_A (int& y, double& a0, double& a1, double& a2, int& m, long double& result, double adjj);
void R0_E1(int& x, int& y, int& k, int& m, double& a0, long double& result);
long double log_R0_E1(int& x, int& y, int& k, int& m, double& a0);
long double	log_R0_E2(int& x, double& a0, int& k);
long double	log_R0_E3(int& y, double& a0, int& m);
void R1_E1(int& k, double& a1, long double& result);
long double log_R1_E1(int& k, double& a1);

void dBvZINB_Expt(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                  double &b1, double &b2, double &p1, double &p2, double &p3, double &p4,
                  Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info, int se);
void dBvZINB_Expt_vec(const Rcpp::IntegerVector &xvec, const Rcpp::IntegerVector &yvec, 
                      const Rcpp::IntegerVector &freq, 
                      const int &n, double &a0, double &a1, double &a2,
                      double &b1, double &b2, double &p1, double &p2, double &p3, double &p4,
                      Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info, int se);

double inv_digamma(double x, double y);
void inv_digamma_vec(Rcpp::NumericVector& lb, Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
                     Rcpp::NumericVector &idgam);
double hfunc(Rcpp::NumericVector& lb, Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
             Rcpp::NumericVector &idgam);
void opt_lb(Rcpp::NumericVector& lb, Rcpp::NumericVector &expt, Rcpp::NumericVector &a,  
            Rcpp::NumericVector &idgam);

void em(Rcpp::NumericVector& param, const Rcpp::IntegerVector &xvec, const Rcpp::IntegerVector &yvec, 
        const Rcpp::IntegerVector &freq, const int &n, Rcpp::NumericVector &expt, Rcpp::NumericVector &info,
        const int &se, Rcpp::IntegerVector &iter, int &maxiter, double &tol, int showFlag);

#endif