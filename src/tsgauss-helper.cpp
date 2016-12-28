#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double tsgauss_lfunc_helper(
    IntegerVector d, NumericVector w, double cutoff,
    NumericVector mu_u, double sd_u, double sigma)
{
  int n = d.size();
  double out = 0;
  for (int i = 0; i < n; i++)
  {
    out += R::pnorm5(w[i] - cutoff, mu_u[i], sd_u, d[i], 1);
  }
  return out / (double)n;
}



// [[Rcpp::export]]
NumericVector tsgauss_lfunc_each_helper(
    IntegerVector d, NumericVector w, double cutoff,
    NumericVector mu_u, double sd_u, double sigma)
{
  int n = d.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++)
  {
    out[i] = R::pnorm5(w[i] - cutoff, mu_u[i], sd_u, d[i], 1);
  }
  return out;
}



/*** R
*/
