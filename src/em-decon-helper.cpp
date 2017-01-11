#include <vector>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double emdecon_update_sigma_gauss(
    double sigma,
    std::vector<double> w_vec, std::vector<int> d_vec, double cutoff,
    std::vector<double> x, std::vector<double> px)
{
  int nobs = w_vec.size();
  int xsize = x.size();
  double h_vec(nobs);
  for (int i = 0; i < nobs; i++)
  {
    std::vector<double> px_pu(xsize);
    for (int j = 0; j < xsize; j++)
      px_pu[j] = px[j] * R::dnorm4(w_vec[i] - x[j], 0, sigma, 0);

    // integrate px_pu for c -> inf or -inf -> c
    // trapezoid approach with fixed point data
    double I1 = 0;
    double I0 = 0;
    for (int j = 0; j < xsize-1; j++)
    {
      if (x[j] <= cutoff && x[j+1] <= cutoff) {
        I0 += (px_pu[j] + px_pu[j+1]);
      } else if (x[j] > cutoff && x[j+1] > cutoff) {
        I1 += (px_pu[j] + px_pu[j+1]);
      } else if (x[j] <= cutoff && x[j+1] > cutoff) {
        double p_mid = px_pu[j]*(x[j+1]-cutoff)/(x[j+1]-x[j]) +
          px_pu[j+1]*(cutoff-x[j])/(x[j+1]-x[j]);
        I0 += (px_pu[j] + p_mid);
        I1 += (px_pu[j+1] + p_mid);
      } else {
        stop("x must be increasing");
      }
    }

  }

  return sigma;
}



/*** R

*/
