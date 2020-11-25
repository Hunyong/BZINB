#ifndef BZINB_H
#define BZINB_H

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <string>
#include <iostream>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
// 
// double EPSILON1;
// double EPSILON2;

void logLin(Rcpp::NumericMatrix &mat, Rcpp::NumericVector &vec, int &n, int &p, int sign, 
            Rcpp::NumericVector &result);
void updateAE(Rcpp::NumericMatrix &ZZ, Rcpp::NumericMatrix &expt, Rcpp::NumericVector &exptSum, Rcpp::NumericVector &b1, 
              Rcpp::NumericVector &eta1, Rcpp::NumericVector &alpha, Rcpp::NumericVector &AEnew, 
              arma::mat &V_eta1, Rcpp::NumericVector &Zbar, arma::mat &score, 
              int &n, int &p, double &error);
void updateEpsilon(Rcpp::NumericMatrix &ZZ, Rcpp::NumericMatrix &expt, Rcpp::NumericVector &b21, 
                   Rcpp::NumericVector &epsilon, Rcpp::NumericVector &epsNew, 
                   arma::mat &V_eps, Rcpp::NumericVector &Zbar, arma::mat &score, 
                   int &n, int &p, double &error);
void updateGamma(Rcpp::NumericMatrix &WW, Rcpp::NumericMatrix &expt,
                 Rcpp::NumericVector &p1, Rcpp::NumericVector &p2, Rcpp::NumericVector &p3, 
                 Rcpp::NumericVector &gamma1, Rcpp::NumericVector &gamma2, Rcpp::NumericVector &gamma3,
                 Rcpp::NumericVector &gammaNew, 
                 arma::mat &V, arma::mat &score, 
                 int &n, int &p, double &error, Rcpp::IntegerVector &zz);
void l1(int& x, int& y, double& a0, double& a1, double& a2, int &k, int& m, 
        double& result, double adjj = 0);
void l1_c (double& t1, double& t2, int& k, int& m, double& result, double adjj);
void l1_AC (double& t1, double& t2, int& x, int& y, double& a0, double& a1, 
            double& a2, int& k, int& m, double& result, double adjj = 0);
void l2_A (int& x, double& a0, double& a1, double& a2, int& k, 
           double& result, double adjj);
void l3_A (int& y, double& a0, double& a1, double& a2, int& m, 
           double& result, double adjj);
void R0_E1(int& x, int& y, int& k, int& m, double& a0, double& result);
double log_R0_E1(int& x, int& y, int& k, int& m, double& a0);
double	log_R0_E2(int& x, double& a0, int& k);
double	log_R0_E3(int& y, double& a0, int& m);
void R1_E1(int& k, double& a1, double& result);
double log_R1_E1(int& k, double& a1);

void dBvZINB_Expt(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                  double &b1, double &b2, double &p1, double &p2,
                  double &p3, double &p4,
                  Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info,
                  int se, int zi, int infoReg);
void dBvZINB_Expt_vec(Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
                      Rcpp::IntegerVector &freq, 
                      int &n, double &a0, double &a1, double &a2,
                      double &b1, double &b2, double &p1, double &p2, 
                      double &p3, double &p4,
                      Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info,
                      int se, int zi);
void dBvZINB_Expt_mat(Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec,
                      Rcpp::NumericMatrix &ZZ, Rcpp::NumericMatrix &WW,
                      int &n, int &pZ, int &pW,
                      Rcpp::NumericVector &alpha, Rcpp::NumericVector &b1, Rcpp::NumericVector &b2, 
                      Rcpp::NumericVector &p1, Rcpp::NumericVector &p2, 
                      Rcpp::NumericVector &p3, Rcpp::NumericVector &p4,
                      Rcpp::NumericMatrix &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &s_i_abp,
                      Rcpp::NumericVector &info, int se, int zi, Rcpp::NumericVector &expt_i);
void dBvZINB_Expt_direct(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                         double &b1, double &b2, double &p1, double &p2, 
                         double &p3, double &p4,
                         Rcpp::NumericVector &l_sum, Rcpp::NumericVector &s_i, 
                         Rcpp::NumericVector &info);
void dBvZINB_Expt_direct_vec(Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
                             Rcpp::IntegerVector &freq, 
                             int &n, double &a0, double &a1, double &a2,
                             double &b1, double &b2, double &p1, double &p2, 
                             double &p3, double &p4,
                             Rcpp::NumericVector &l_sum, Rcpp::NumericVector &s_i, 
                             Rcpp::NumericVector &info);

double inv_digamma(double x, double y);
void inv_digamma_vec(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
                     double idgam[3]);
double hfunc(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
                  double idgam[3]);
void opt_lb(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
            double idgam[3]);

Rcpp::List em(Rcpp::NumericVector& param2, Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
              Rcpp::IntegerVector &freq, int &n, int &se, int &maxiter, double &tol, int showFlag, 
              int zi);

Rcpp::List emReg(Rcpp::NumericVector& param2, Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
                 Rcpp::NumericMatrix& ZZ, Rcpp::NumericMatrix& WW,
                 int &pZ, int &pW,
                 int &n, int &se, int &maxiter, double &tol, int showFlag,
                 Rcpp::IntegerVector &zi);

#endif