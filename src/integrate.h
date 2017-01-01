#ifndef INTEGRATEHEADERDEF
#define INTEGRATEHEADERDEF

#include <Rcpp.h>
#include <functional>
#include <cmath>
#include <string>
using namespace Rcpp;


// numerical integration routine for univariate real function


// function 'Integrate' is the wrapper for individual method
double Integrate(std::function<double(double)> &func,
                 double a, double b,
                 std::string method, double tol);


// adaptive integration by trapezoid method
double IntegrateTrapezoid(std::function<double(double)> &func,
                          double a, double b, double tol);
double IntegrateTrapezoid(std::function<double(double)> &func,
                          double a, double b, double f_a, double f_b,
                          double tol);


// adoptive integration by simpson's 1/3 method
double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double tol);
double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double mid,
                        double f_a, double f_b, double f_mid,
                        double tol);


// adoptive integration by simpson's 3/8 method
double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double tol);
double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double m1, double m2,
                         double f_a, double f_b, double f_m1, double f_m2,
                         double tol);


// adoptive integration by Romberg method
double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double tol);
double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double mid,
                        double f_a, double f_b, double f_mid,
                        double tol);

#endif