#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <string>
#include "integral.h"

using namespace Rcpp;


class EmGaussLapModel
{

private:
  // data
  std::vector<int> d_vec;
  std::vector<double> w_vec;
  int nobs;
  double cutoff;

  // control parameters
  double integ_tol;  // accuracy for integration error
  int integ_depth;   // maximum recursion depth for integration
  std::string integ_method; // integration method
  double tol;        // threshold tolerance for value improvement
  int maxit;

  // model parameters
  double sigma;
  double mu_x;
  double sd_x;
  double sd_w;

  // and asymptotic variance and standard errors
  NumericMatrix avar;
  NumericVector stderr;

  // value and weights, and convergence indicator
  double cur_value;
  std::vector<double> weights;
  int convergence;

  // density functions
  double fx(double x)
  { return R::dnorm4(x, mu_x, sd_x, 0); }
  double fu(double x, int i)
  { return exp(-sqrt(2) / sigma * fabs(w_vec[i] - x)) / sigma / sqrt(2.0); }



  // returns increment
  double UpdateValueAndWeights()
  {
    double new_value = 0;
    for (int i = 0; i < nobs; i++)
    {
      std::function<double(double)> func = [this, i] (double x) -> double {
        return fx(x) * fu(x, i); };

      double lower;
      double upper;
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      weights[i] = Integrate(func, lower, upper,
                             integ_method, integ_tol, integ_depth);
      new_value += log(weights[i]);
    }
    new_value /= double (nobs);
    double increment = new_value - cur_value;
    cur_value = new_value;
    return increment;
  }

  void UpdateParameters()
  {
    // update parameters for u
    std::function<double(double)> func;
    double new_sigma = 0;
    for (int i = 0; i < nobs; i++)
    {
      func = [this,i] (double x) -> double {
        return fx(x) * fu(x, i) / weights[i] * fabs(w_vec[i] - x); };
      double lower;
      double upper;
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      new_sigma += Integrate(func, lower, upper,
                             integ_method, integ_tol, integ_depth);
    }
    new_sigma *= (sqrt(2.0) / (double)nobs);

    // update paramters for x
    double new_sdx = 0;
    for (int i = 0; i < nobs; i++)
    {
      func = [this,i] (double x) -> double {
        return fx(x) * fu(x, i) / weights[i] * pow(x - mu_x, 2); };
      double lower;
      double upper;
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      new_sdx += Integrate(func, lower, upper,
                             integ_method, integ_tol, integ_depth);
    }
    new_sdx /= (double)nobs;
    new_sdx = sqrt(new_sdx);

    sigma = new_sigma;
    sd_x = new_sdx;
  }


  void UpdateAvarAndSe()
  {
    // implements Murphy and Topel (1985), section 5.1
    arma::mat R1(1, 1);
    arma::mat R2(2, 2);
    arma::mat R3(1, 2);
    arma::mat R4(1, 2);
    R1(0) = 1.0 / sd_w;

    // compute Jacobians, namely score vectors
    arma::mat J11(nobs, 1);  // diff L1 on theta1
    arma::mat J21(nobs, 1);  // diff L2 on theta1
    arma::mat J22(nobs, 2);  // diff L2 on theta2
    for (int i = 0; i < nobs; i++)
    {
      // have analytic solution for J11, derived from gaussian pdf
      J11(i, 0) = (w_vec[i] - mu_x)/pow(sd_w, 2);

      // J21 and J22 requires numerical integration
      std::function<double(double)> func;
      double lower;
      double upper;

      // computing L2 on mu_x
      func = [this,i] (double x) -> double {
        return (x - mu_x) / pow(sd_x, 2) * fx(x) * fu(x, i); };
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      J21(i, 0) = Integrate(func, lower, upper,
          integ_method, integ_tol, integ_depth) / weights[i];

      // computing L2 on sigma
      func = [this,i] (double x) -> double {
        return (-1.0/sigma + sqrt(2)*fabs(w_vec[i]-x)/pow(sigma, 2)) *
          fx(x) * fu(x, i); };
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      J22(i, 0) = Integrate(func, lower, upper,
          integ_method, integ_tol, integ_depth) / weights[i];

      // computing L2 on sd_x
      func = [this,i] (double x) -> double {
        return (-1.0/sd_x + pow(x - mu_x, 2)/pow(sd_x, 3)) *
          fx(x) * fu(x, i); };
      if (d_vec[i] == 1) {
        lower = cutoff;
        upper = INFINITY;
      } else {
        lower = -INFINITY;
        upper = cutoff;
      }
      J22(i, 1) = Integrate(func, lower, upper,
          integ_method, integ_tol, integ_depth) / weights[i];
    }
    R2 = J22.t() * J22 / nobs;
    R3 = J21.t() * J22 / nobs;
    R4 = J11.t() * J22 / nobs;


    arma::mat o11 = R1.i();
    arma::mat tmp = R2.i();
    arma::mat o12 = o11 * (R4-R3) * tmp;
    arma::mat o22 = tmp + tmp *
       (R3.t() * o11 * R3 - R4.t() * o11 * R3 - R3.t() * o11 * R4) * tmp;
    arma::mat o21 = o12.t();

    arma::mat out = arma::join_vert(arma::join_horiz(o11, o12),
                                    arma::join_horiz(o21, o22));

    // need to arrange matrix
    std::vector<int> ind = {1, 0, 2};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        avar(i,j) = out(ind[i], ind[j]);

    // update standard errors
    for (int i = 0; i < 3; i++)
      stderr[i] = sqrt(avar(i,i)/nobs);

  }


public:

  EmGaussLapModel(
    const std::vector<int> &d_vec_, const std::vector<double> &w_vec_,
    double cutoff_, double tol_, int maxit_,
    std::string integ_method_, double integ_tol_, int integ_depth_)
  {
    d_vec = d_vec_;
    w_vec = w_vec_;
    nobs = d_vec.size();
    cutoff = cutoff_;
    tol = tol_;
    maxit = maxit_,
    integ_tol = integ_tol_;
    integ_depth = integ_depth_;
    integ_method = integ_method_;

    // initialize parameters
    // set mu_x = mean(w)
    // initialize sd_x = sigma = sqrt(var(w)/2)
    mu_x = 0.0;
    double w2 = 0.0;
    for (int i = 0; i < nobs; i++)
    {
      mu_x += w_vec[i];
      w2 += w_vec[i]*w_vec[i];
    }
    mu_x /= (double) nobs;
    double s2 = w2 / (double) nobs - mu_x*mu_x; // (sigma_w)^2
    sd_x = sqrt(s2 * 0.75);
    sigma = sqrt(s2 * 0.25);
    sd_w = sqrt(s2);

    // initialize avar and se
    avar = NumericMatrix(3, 3);
    rownames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    colnames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    stderr = NumericVector::create(
      Named("sigma") = 0, Named("mu_x") = 0, Named("sd_x") = 0);


    // initialize value and weights
    cur_value = -INFINITY;
    weights.resize(nobs, 1);
    UpdateValueAndWeights();


    convergence = -1;
  }

  void Estimate(bool verbose)
  {
    convergence = 1;
    for (int i = 0; i < maxit; i++)
    {
      if (verbose) {
        Rcout.precision(3);
        Rcout << "iter " << i+1 <<
          " sigma = " << sigma << ", sd_x = " << sd_x <<
          " value = " << cur_value << "\n";
      }
      UpdateParameters();
      double increment = UpdateValueAndWeights();
      if (fabs(increment) < tol*(fabs(cur_value) + tol)) {
        convergence = 0;
        if (verbose) Rcout << "CONVERGED!\n";
        break;
      }
    }

    UpdateAvarAndSe();
  }

  List CompileOutput()
  {
    NumericVector estimate = NumericVector::create(
      Named("sigma") = sigma, Named("mu_x") = mu_x, Named("sd_x") = sd_x
    );

    List out = List::create(
      Named("estimate") = estimate,
      Named("stderr") = stderr,
      Named("avar") = avar,
      Named("nobs") = nobs,
      Named("convergence") = convergence
    );
    return out;

  }

  void Summary()
  {
    Rcout << "Parameters:\n" <<
      "  sigma = " << sigma << "\n" <<
      "  mu_x  = " << mu_x << "\n" <<
      "  sd_x  = " << sd_x << "\n" <<
      "Value = " << cur_value << "\n\n";
  }

};


// [[Rcpp::export]]
List em_gauss_lap_helper(
    std::vector<int> d_vec, std::vector<double> w_vec, double cutoff,
    double tol, int maxit,
    std::string integ_method, double integ_tol, int integ_depth,
    bool verbose)
{
  EmGaussLapModel model(d_vec, w_vec, cutoff, tol, maxit,
                        integ_method, integ_tol, integ_depth);
  model.Estimate(verbose);
  return model.CompileOutput();
}


/*** R
library(rddsigma)
dat <- gen_data(500, 0.3, 1)
a <- Sys.time()
rddsigma:::em_gauss_lap_helper(dat$d, dat$w, 1,
                               1e-6, 1000, "romberg", 1e-6, 100)
Sys.time() - a
a <- Sys.time()
em_gauss_lap(dat$d, dat$w, 1, reltol = 1e-6,
             integrate_options = list(rel.tol = 1e-6), quiet=TRUE)
Sys.time() - a
*/
