#include <Rcpp.h>
#include "integrate.h"
using namespace Rcpp;


double Integrate(std::function<double(double)> &func,
                 double a, double b,
                 std::string method, double tol, int max_depth)
{
  // handles infinity in a and/or b by change of variable approach
  if (isinf(a) && isinf(b)) {
    // if loewr = +inf or upper = -inf --> integral is 0
    if (a > 0 || b < 0) return 0.0;

    // split into two, (-inf, 0) and (0, +inf)
    return Integrate(func, a, 0.0, method, tol, max_depth) +
      Integrate(func, 0.0, b, method, tol, max_depth);
  } else if (isinf(a)) {
    if (a > 0) return 0.0;  // lower is +inf --> integral is 0

    // this is an integral over (-inf, b). to do this,
    // change the variable to t = -x and integrate over (-b, +inf)
    std::function<double(double)> new_func = [&func](double z) -> double {
      return func(-z); };
    return Integrate(new_func, -b, INFINITY, method, tol, max_depth);
  } else if (isinf(b)) {
    if (b < 0) return 0.0;  // upper if -inf --> integral is 0

    // this is an integral over (a, +inf).
    // we apply a change of variable:
    // z = (x-a) / (1+x-a) or x = z/(1-z) + a
    // then,
    // I = integral(0, 1) 1/(1-z)^2 f(z/(1-z) + a) dz
    std::function<double(double)> new_func = [&func, a](double z) -> double {
      return func(z/(1-z) + a) * pow(1-z, -2.0); };

    if (fabs(func(b)) > 0.0) stop("function does not converge to zero");

    return Integrate(new_func, 0.0, 1.0, new_func(0.0), 0.0,
                     method, tol, max_depth);
    //return Integrate(new_func, 0.0, 1.0, method, tol, max_depth);
  }

  return Integrate(func, a, b, func(a), func(b), method, tol, max_depth);
}


double Integrate(std::function<double(double)> &func,
                 double a, double b, double f_a, double f_b,
                 std::string method, double tol, int max_depth)
{
  if (method == "trapezoid") {
    return IntegrateTrapezoid(func, a, b, f_a, f_b, tol, max_depth);
  } else if (method == "simpson") {
    return IntegrateSimpson(func, a, b, f_a, f_b, tol, max_depth);
  } else if (method == "simpson2") {
    return IntegrateSimpson2(func, a, b, f_a, f_b, tol, max_depth);
  } else if (method == "romberg") {
    return IntegrateRomberg(func, a, b, f_a, f_b, tol, max_depth);
  }
  stop("unknown integration method");
  return 0;
}



// double IntegrateTrapezoid(std::function<double(double)> &func,
//                           double a, double b, double f_a, double f_b,
//                           double tol, int max_depth)
// {
//   return IntegrateTrapezoid(func, a, b, f_a, f_b, tol, max_depth);
// }

double IntegrateTrapezoid(std::function<double(double)> &func,
                          double a, double b, double f_a, double f_b,
                          double tol, int depth_remained)
{
  double mid = 0.5*(a + b);
  double f_mid = func(mid);

  // stop if nan or inf is observed
  if (isnan(f_a) || isnan(f_b) || isnan(f_mid))
    stop("point value is nan");
  if (isinf(f_a) || isinf(f_b) || isinf(f_mid))
    stop("point value is inf or -inf");

  double I1 = 0.5*(f_a + f_b)*(b-a);
  double I2 = 0.5*(f_a + 2*f_mid + f_b)*0.5*(b-a);

  if (fabs(I2 - I1) < tol*(fabs(I1) + tol)) return I2;

  depth_remained--;
  if (depth_remained < 0) stop("maximum recursion depth exceeded");

  return IntegrateTrapezoid(func, a, mid, f_a, f_mid, tol, depth_remained) +
    IntegrateTrapezoid(func, mid, b, f_mid, f_b, tol, depth_remained);
}



double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double f_a, double f_b,
                        double tol, int max_depth)
{
  double mid = 0.5*(a+b);
  return IntegrateSimpson(func, a, b, mid, f_a, f_b, func(mid),
                          tol, max_depth);
}

double IntegrateSimpson(std::function<double(double)> &func,
                        double a, double b, double mid,
                        double f_a, double f_b, double f_mid,
                        double tol, int depth_remained)
{
  double m1 = 0.5*(a+mid);
  double m2 = 0.5*(mid+b);
  double f_m1 = func(m1);
  double f_m2 = func(m2);

  // stop if nan or inf is observed
  if (isnan(f_a) || isnan(f_b) || isnan(f_m1) || isnan(f_m2))
    stop("point value is nan");
  if (isinf(f_a) || isinf(f_b) || isinf(f_m1) || isinf(f_m1))
    stop("point value is inf or -inf");

  double I1 = (f_a + 4*f_mid + f_b) / 3.0 * 0.5*(b-a);
  double I2 = (f_a + 4.0*f_m1 + 2.0*f_mid + 4.0*f_m2 + f_b) / 3.0 * 0.25*(b-a);

  if (fabs(I2-I1) < tol*(fabs(I1) + tol)) return I2;

  depth_remained--;
  if (depth_remained < 0) stop("maximum recursion depth exceeded");
  return IntegrateSimpson(func, a, mid, m1, f_a, f_mid, f_m1, tol, depth_remained) +
    IntegrateSimpson(func, mid, b, m2, f_mid, f_b, f_m2, tol, depth_remained);
}


double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double f_a, double f_b,
                         double tol, int max_depth)
{
  double m1 = a*2.0/3.0 + b/3.0;
  double m2 = a/3.0 + b*2.0/3.0;
  return IntegrateSimpson2(func, a, b, m1, m2,
                           f_a, f_b, func(m1), func(m2),
                           tol, max_depth);
}

double IntegrateSimpson2(std::function<double(double)> &func,
                         double a, double b, double m1, double m2,
                         double f_a, double f_b, double f_m1, double f_m2,
                         double tol, int depth_remained)
{
  double m0 = 0.5*(a+m1);
  double mid = 0.5*(m1+m2);
  double m3 = 0.5*(m2+b);
  double f_m0 = func(m0);
  double f_mid = func(mid);
  double f_m3 = func(m3);

  // stop if nan or inf is observed
  if (isnan(f_a) || isnan(f_b) || isnan(f_mid) || isnan(f_m0) || isnan(f_m3))
    stop("point value is nan");
  if (isinf(f_a) || isinf(f_b) || isinf(f_mid) || isinf(f_m0) || isinf(f_m3))
    stop("point value is inf or -inf");

  double I1 = 3.0/8.0 * (b-a)/3.0 * (f_a + 3.0*f_m1 + 3.0*f_m1 + f_b);
  double I2 = 3.0/8.0 * (b-a)/6.0 *
    (f_a + 3.0*f_m0 + 3.0*f_m1 + 2.0*f_mid + 3.0*f_m2 + 3.0*f_m3 + f_b);

  if (fabs(I2-I1) < tol*(fabs(I1) + tol)) return I2;

  depth_remained--;
  if (depth_remained < 0) stop("maximum recursion depth exceeded");
  return IntegrateSimpson2(func, a, mid, m0, m1, f_a, f_mid, f_m0, f_m1,
                           tol, depth_remained) +
    IntegrateSimpson2(func, mid, b, m2, m3, f_mid, f_b, f_m2, f_m3,
                      tol, depth_remained);
}




double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double f_a, double f_b,
                        double tol, int max_depth)
{
  double mid = 0.5*(a+b);
  return IntegrateRomberg(func, a, b, mid, f_a, f_b, func(mid),
                          tol, max_depth);
}

double IntegrateRomberg(std::function<double(double)> &func,
                        double a, double b, double mid,
                        double f_a, double f_b, double f_mid,
                        double tol, int depth_remained)
{
  double m1 = 0.5*(a+mid);
  double m2 = 0.5*(mid+b);
  double f_m1 = func(m1);
  double f_m2 = func(m2);

  // stop if nan or inf is observed
  if (isnan(f_a) || isnan(f_b) || isnan(f_mid) || isnan(f_m1) || isnan(f_m2))
    stop("point value is nan");
  if (isinf(f_a) || isinf(f_b) || isinf(f_mid) || isinf(f_m1) || isnan(f_m2))
    stop("point value is inf or -inf");

  double I1 = 0.5 * (b-a) * (f_a + f_b);
  double I2 = 0.5 * 0.5*(b-a) * (f_a + 2.0*f_mid + f_b);
  double I3 = 0.5 * 0.25*(b-a) * (f_a + 2.0*f_m1 + 2.0*f_mid + 2.0*f_m2 + f_b);

  double I4 = I2 + (I2 - I1)/3.0;
  double I5 = I3 + (I3 - I2)/3.0;

  if (fabs(I5-I4) < tol*(fabs(I4) + tol)) return I5;

  depth_remained--;
  if (depth_remained < 0) stop("maximum recursion depth exceeded");
  return IntegrateRomberg(func, a, mid, m1, f_a, f_mid, f_m1,
                          tol, depth_remained) +
    IntegrateRomberg(func, mid, b, m2, f_mid, f_b, f_m2,
                     tol, depth_remained);
}



// [[Rcpp::export]]
void integrate_test()
{
  Rcout.precision(10);

  double tol = 1e-5;
  int max_depth = 100;
  std::function<double(double)> f1 = [](double x) -> double { return 1/x; };
  Rcout << "integrate 1/x over (3.1, 3.9)\n" <<
    "  true value = " << log(3.9) - log(3.1) << "\n" <<
    "  trapezoid  = " << Integrate(f1, 3.1, 3.9, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f1, 3.1, 3.9, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f1, 3.1, 3.9, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f1, 3.1, 3.9, "romberg", tol, max_depth) << "\n\n";

  std::function<double(double)> f2 = [](double x) -> double {
    return exp(-0.5*x*x)/sqrt(2*3.141592653589793); };
  Rcout << "integrate standard normal pdf over (-1.96, 1.96)\n" <<
    "  true value = " <<
    R::pnorm5(1.96, 0.0, 1.0, 1, 0) - R::pnorm5(-1.96, 0.0, 1.0, 1, 0) << "\n" <<
    "  trapezoid  = " << Integrate(f2, -1.96, 1.96, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, -1.96, 1.96, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, -1.96, 1.96, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, -1.96, 1.96, "romberg", tol, max_depth) << "\n\n";

  Rcout << "integrate standard normal pdf over (0, +inf)\n" <<
    "  true value = 0.5\n" <<
    "  trapezoid  = " << Integrate(f2, 0, INFINITY, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, 0, INFINITY, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, 0, INFINITY, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, 0, INFINITY, "romberg", tol, max_depth) << "\n\n";

  Rcout << "integrate standard normal pdf over (-inf, 0)\n" <<
    "  true value = 0.5\n" <<
    "  trapezoid  = " << Integrate(f2, -INFINITY, 0, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, -INFINITY, 0, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, -INFINITY, 0, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, -INFINITY, 0, "romberg", tol, max_depth) << "\n\n";

  Rcout << "integrate standard normal pdf over (-inf, +inf)\n" <<
    "  true value = 1.0\n" <<
    "  trapezoid  = " << Integrate(f2, -INFINITY, INFINITY, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, -INFINITY, INFINITY, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, -INFINITY, INFINITY, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, -INFINITY, INFINITY, "romberg", tol, max_depth) << "\n\n";

  Rcout << "integrate standard normal pdf over (-1.64, +inf)\n" <<
    "  true value = " << R::pnorm5(-1.64, 0.0, 1.0, 0, 0) << "\n"
    "  trapezoid  = " << Integrate(f2, -1.64, INFINITY, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, -1.64, INFINITY, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, -1.64, INFINITY, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, -1.64, INFINITY, "romberg", tol, max_depth) << "\n\n";

  Rcout << "integrate standard normal pdf over (-inf, 1.28)\n" <<
    "  true value = " << R::pnorm5(1.28, 0.0, 1.0, 1, 0) << "\n"
    "  trapezoid  = " << Integrate(f2, -INFINITY, 1.28, "trapezoid", tol, max_depth) << "\n" <<
    "  simpson    = " << Integrate(f2, -INFINITY, 1.28, "simpson", tol, max_depth) << "\n" <<
    "  simpson2   = " << Integrate(f2, -INFINITY, 1.28, "simpson2", tol, max_depth) << "\n" <<
    "  romberg    = " << Integrate(f2, -INFINITY, 1.28, "romberg", tol, max_depth) << "\n\n";
}

// [[Rcpp::export]]
void integrate_test2(std::string method)
{
  // compute (-inf, inf) integral for standard normal pdf (the ans is 1.0)
  // this is for comparing the speed
  std::function<double(double)> f = [](double x) -> double {
    return exp(-0.5*x*x)/sqrt(2*3.141592653589793); };
  Integrate(f, -INFINITY, INFINITY, method, 1e-5, 100);
}

// [[Rcpp::export]]
void integrate_test3()
{
  // compute divergent integral
  // this would create max recursion error
  std::function<double(double)> f = [](double x) -> double {
    return 1/x; };
  Rcout.precision(10);
  Rcout << "integrate 1/x over (1, 1e+10)\n" <<
    "  trapezoid: " << Integrate(f, 1, 1e+10, "trapezoid", 1e-5, 100) << "\n" <<
    "  simpson:   " << Integrate(f, 1, 1e+10, "simpson", 1e-5, 100) << "\n" <<
    "  simpson2:  " << Integrate(f, 1, 1e+10, "simpson2", 1e-5, 100) << "\n" <<
    "  romberg:   " << Integrate(f, 1, 1e+10, "romberg", 1e-5, 100) << "\n";
  Rcout << "integrate 1/x over (1, 1e+20)\n" <<
    "  trapezoid: " << Integrate(f, 1, 1e+20, "trapezoid", 1e-5, 100) << "\n" <<
    "  simpson:   " << Integrate(f, 1, 1e+20, "simpson", 1e-5, 100) << "\n" <<
    "  simpson2:  " << Integrate(f, 1, 1e+20, "simpson2", 1e-5, 100) << "\n" <<
    "  romberg:   " << Integrate(f, 1, 1e+20, "romberg", 1e-5, 100) << "\n";
  Rcout << "integrate 1/x over (1, +inf)\n" <<
    "  trapezoid: " << Integrate(f, 1, INFINITY, "trapezoid", 1e-5, 100) << "\n" <<
    "  simpson:   " << Integrate(f, 1, INFINITY, "simpson", 1e-5, 100) << "\n" <<
    "  simpson2:  " << Integrate(f, 1, INFINITY, "simpson2", 1e-5, 100) << "\n" <<
    "  romberg:   " << Integrate(f, 1, INFINITY, "romberg", 1e-5, 100) << "\n";
}



/*** R
rddsigma:::integrate_test()
library(microbenchmark)
library(ggplot2)
m <- microbenchmark(
  trapezoid = rddsigma:::integrate_test2("trapezoid"),
  simpson   = rddsigma:::integrate_test2("simpson"),
  simpson2  = rddsigma:::integrate_test2("simpson2"),
  romberg   = rddsigma:::integrate_test2("romberg")
)
summary(m)
ggplot(as.data.frame(m), aes(expr, time)) + geom_boxplot() + scale_y_log10()

rddsigma:::integrate_test3()
*/
