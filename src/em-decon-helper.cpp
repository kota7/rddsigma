#include <vector>
#include <Rcpp.h>
using namespace Rcpp;



double integration_helper(
    std::vector<double> &x, std::vector<double> &fx,
    double c, bool lt)
{
  // This function computes integral from fixed point values
  //
  // Args:
  //   x  : vector of evaluation point, assumed to be increasing
  //   fx : vector of function values
  //   c  : integration limit
  //   lt : compute lower tail or upper tail
  //
  // Returns:
  //   if lt is true, then
  //     integral_{-inf, c} f(x) dx
  //   if lt is false, then
  //     integral_{c, +inf} f(x) dx
  //
  // Note:
  //   implements trapezoid approach
  if (x.size() != fx.size()) stop("x and fx must have the same length");
  if (x.size() < 2) stop("size of x is too small");

  int xsize = x.size();
  double I = 0;
  if (lt) {
    for (int j = 0; j < xsize-1; j++)
    {
      if (x[j] <= c && x[j+1] <= c) {
        I += (fx[j] + fx[j+1]);
      } else if (x[j] <= c && x[j+1] > c) {
        double p_mid = fx[j]*(x[j+1]-c)/(x[j+1]-x[j]) +
          fx[j+1]*(c-x[j])/(x[j+1]-x[j]);
        I += (fx[j] + p_mid);
        break;
      }
    }
  } else {
    for (int j = xsize-2; j >= 0; j--)
    {
      if (x[j] > c && x[j+1] > c) {
        I += (fx[j] + fx[j+1]);
      } else if (x[j] <= c && x[j+1] > c) {
        double p_mid = fx[j]*(x[j+1]-c)/(x[j+1]-x[j]) +
          fx[j+1]*(c-x[j])/(x[j+1]-x[j]);
        I += (fx[j+1] + p_mid);
        break;
      }
    }
  }
  return I;
}


// [[Rcpp::export]]
double emdecon_update_sigma_gauss(
    double sigma,
    std::vector<int> d_vec, std::vector<double> w_vec, double cutoff,
    std::vector<double> x, std::vector<double> px)
{
  // input validation
  // d must be either 1 or 0
  for (std::vector<int>::iterator it = d_vec.begin(); it != d_vec.end(); ++it)
  {
    if (*it != 1 && *it != 0) stop("d must be vector of 1 or 0");
  }

  // x must be increasing
  for (size_t i = 0; i < x.size()-1; i++)
  {
    if (x[i] >= x[i+1]) stop("x must be strictly increasing");
  }

  // x and px must be the same length
  if (x.size() != px.size()) stop("x and px must be of the same size");

  // w and d must be the same length
  if (d_vec.size() != w_vec.size()) stop("d and w must be of the same size");



  int nobs = w_vec.size();
  int xsize = x.size();
  double new_sigma = 0;
  for (int i = 0; i < nobs; i++)
  {
    std::vector<double> px_pu(xsize);
    for (int j = 0; j < xsize; j++)
    {
      // this part varies by the distribution of u
      px_pu[j] = px[j] * R::dnorm4(w_vec[i] - x[j], 0, sigma, 0);
    }

    // computation of h_i

    // first, compute the denominator
    // integrate px_pu for c -> inf or -inf -> c
    double I = integration_helper(x, px_pu, cutoff, d_vec[i] == 0);
    // raise error if the integral,
    // i.e. the prob of assignment is not positive
    if (I <= 0) stop("assignment probability is not positive");


    // compute the contribution to new_sigma from observation i
    // this part varies by the distribution of u
    std::vector<double> fx(xsize);
    for (int j = 0; j < xsize; j++)
    {
      // this part varies by the distribution of u
      fx[j] = px_pu[j] / I * std::pow(w_vec[i] - x[j], 2);
    }
    double cont = integration_helper(x, fx, cutoff,  d_vec[i] == 0);
    new_sigma += cont;
  }
  new_sigma = std::sqrt(new_sigma / (double)nobs);

  return new_sigma;
}





// [[Rcpp::export]]
double emdecon_update_sigma_lap(
    double sigma,
    std::vector<int> d_vec, std::vector<double> w_vec, double cutoff,
    std::vector<double> x, std::vector<double> px)
{
  // input validation
  // d must be either 1 or 0
  for (std::vector<int>::iterator it = d_vec.begin(); it != d_vec.end(); ++it)
  {
    if (*it != 1 && *it != 0) stop("d must be vector of 1 or 0");
  }

  // x must be increasing
  for (size_t i = 0; i < x.size()-1; i++)
  {
    if (x[i] >= x[i+1]) stop("x must be strictly increasing");
  }

  // x and px must be the same length
  if (x.size() != px.size()) stop("x and px must be of the same size");

  // w and d must be the same length
  if (d_vec.size() != w_vec.size()) stop("d and w must be of the same size");



  int nobs = w_vec.size();
  int xsize = x.size();
  double new_sigma = 0;
  for (int i = 0; i < nobs; i++)
  {
    std::vector<double> px_pu(xsize);
    for (int j = 0; j < xsize; j++)
    {
      // this part varies by the distribution of u
      px_pu[j] = px[j] *
        exp(-std::sqrt(2.0) / sigma * std::fabs(w_vec[i] - x[j])) / sigma / std::sqrt(2.0);
    }

    // computation of h_i

    // first, compute the denominator
    // integrate px_pu for c -> inf or -inf -> c
    double I = integration_helper(x, px_pu, cutoff, d_vec[i] == 0);
    // raise error if the integral,
    // i.e. the prob of assignment is not positive
    if (I <= 0) stop("assignment probability is not positive");


    // compute the contribution to new_sigma from observation i
    // this part varies by the distribution of u
    std::vector<double> fx(xsize);
    for (int j = 0; j < xsize; j++)
    {
      // this part varies by the distribution of u
      fx[j] = px_pu[j] / I * std::fabs(w_vec[i] - x[j]);
    }
    double cont = integration_helper(x, fx, cutoff, d_vec[i] == 0);
    new_sigma += cont;
  }
  new_sigma = std::sqrt(2.0) * new_sigma / (double)nobs;


  return new_sigma;
}


/*** R

*/
