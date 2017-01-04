#ifndef INTEGRATEHEADERDEF
#define INTEGRATEHEADERDEF

#include <Rcpp.h>
#include <functional>
#include <cmath>
#include <string>

using namespace Rcpp;


// Numerical integration routine for univariate real function
//
// Implements adaptive procedure as below.
// Integrate(a, b):
//   I1 = estimated integral
//   I2 = estimated integral one step finer
//   if abs(I2-I1)/(abs(I1) + rol) < tol:
//     return I2
//   else:
//     return Integrate(a, (a+b)/2) + Integrate((a+b)/2, b)
//
// It is designed so that I2 in the current call coincides
// with I1 in the recursive call.
//
// 4 methodologies are implemented
//   Trapezoid
//   Simpson 1/3
//   Simpson 3/8
//   Romberg (Trapezoid with two-step extrapolation)
// They differ in the way computing I1 and I2
// Reference:
//   Joe D Hoffman (2001) Numerical Methods for Engineers and Scientists.
//     Marcel Dekker, Inc.
//
// Integrand must return a finite value for each point on the range [a, b]
// Error is raised when inf or nan is evaluated.
//
// Improper integrals (lower and/or upper limit is infinite) are handled by
// change of variable:
// z = (x-a) / (1+x-a) or x = z/(1-z) + a
//   then,
// Integral_(a, inf) f(x) dx = integral_(0, 1) 1/(1-z)^2 f(z/(1-z) + a) dz
//
// If the fuction is applied to compute undefined integrals,
// that is, answer is infinite, it would produce 'maximum recursion error'




// this is the function to be used.
double Integrate(std::function<double(double)> &func,
                 double a, double b,
                 std::string method, double tol, int max_depth);
// this version requires function values at the both end.
// this is designed for computing improper integrals
// for which func(a) or func(b) may not be evaluated properly
double Integrate(std::function<double(double)> &func,
                 double a, double b, double f_a, double f_b,
                 std::string method, double tol, int max_depth);






// adaptive integration by trapezoid method
// note: the first version is deleted since it has the same signature as
// the recursive version
// double IntegrateTrapezoid(std::function<double(double)> &func,
//                           double a, double b, double f_a, double f_b,
//                           double tol, int max_depth);
double IntegrateTrapezoid(std::function<double(double)> &func,
                          double a, double b, double f_a, double f_b,
                          double tol, int depth_remained);


// adoptive integration by simpson's 1/3 method
double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double f_a, double f_b,
                        double tol, int max_depth);
double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double mid,
                        double f_a, double f_b, double f_mid,
                        double tol, int depth_remained);


// adoptive integration by simpson's 3/8 method
double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double f_a, double f_b,
                         double tol, int max_depth);
double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double m1, double m2,
                         double f_a, double f_b, double f_m1, double f_m2,
                         double tol, int depth_remained);


// adoptive integration by Romberg method (trapezoid + 2-step extrapolation)
double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double f_a, double f_b,
                        double tol, int max_depth);
double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double m1, double m2, double m3,
                        double f_a, double f_b,
                        double f_m1, double f_m2, double f_m3,
                        double tol, int depth_remained);

#endif