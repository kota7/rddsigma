#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <string>
#include "integral.h"

using namespace Rcpp;


class EmGaussGaussModel
{
private:
  // data
  std::vector<int> d_vec;
  std::vector<double> w_vec;
  int nobs;
  double cutoff;

  // control parameters
  double integ_tol;   // accuracy for integration error
  double integ_tol0;  // accuracy for integration error specified originally
  int integ_depth;   // maximum recursion depth for integration
  std::string integ_method; // integration method
  double tol;        // threshold tolerance for value improvement
  int maxit;

  // error catcher
  double integ_tol_limit;
  bool accuracy_limit_error;

  // model parameters
  double sigma;
  double mu_x;
  double sd_x;
  double sd_w;

  // and asymptotic variance and standard errors
  NumericMatrix avar;
  NumericVector stderror;

  // value and weights, and convergence indicator
  double cur_value;
  std::vector<double> weights;
  int convergence;
  double last_increment;

  // density functions
  double fx(double x)
  { return R::dnorm4(x, mu_x, sd_x, 0); }
  double fu(double x, int i)
  { return R::dnorm4(w_vec[i] - x, 0, sigma, 0); }



  double UpdateValueAndWeights()
  {
    // update value and weights, then
    // returns increment
    double new_value = 0;
    std::vector<double> new_weights(weights.size());

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
      new_weights[i] = Integrate(func, lower, upper,
                                 integ_method, integ_tol, integ_depth);
      new_value += log(new_weights[i]);
    }
    new_value /= double (nobs);
    double increment = new_value - cur_value;

    // check if the increment is positive
    // otherwise, increase integration accuracy and re-compute
    // but we allow tiny negative increment
    if (increment < -tol*(std::fabs(cur_value) + tol)) {
      integ_tol *= 0.2;
      Rcout << "negative increment: " << increment <<
          " integ_tol is now " << integ_tol << "\n";

      // terminate if the accuracy is larger than the required limit
      if (integ_tol < integ_tol_limit) {
        accuracy_limit_error = true;
        return increment;
      }

      // recompute parameters with new accuracy
      UpdateParameters();
      return UpdateValueAndWeights();
    }

    // update value
    cur_value = new_value;
    for (int i = 0; i < weights.size(); i++) weights[i] = new_weights[i];

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
        return fx(x) * fu(x, i) / weights[i] * std::pow(w_vec[i] - x, 2); };
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
    new_sigma = std::sqrt(new_sigma / (double)nobs);

    // update paramters for x
    double new_sdx = 0;
    for (int i = 0; i < nobs; i++)
    {
      func = [this,i] (double x) -> double {
        return fx(x) * fu(x, i) / weights[i] * std::pow(x - mu_x, 2); };
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
    new_sdx = std::sqrt(new_sdx / (double)nobs);

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
      J11(i, 0) = (w_vec[i] - mu_x)/std::pow(sd_w, 2);

      // J21 and J22 requires numerical integration
      std::function<double(double)> func;
      double lower;
      double upper;

      // computing L2 on mu_x
      func = [this,i] (double x) -> double {
        return (x - mu_x) / std::pow(sd_x, 2) * fx(x) * fu(x, i); };
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
        return (-1.0/sigma + std::pow(w_vec[i]-x, 2)/std::pow(sigma, 3)) *
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
        return (-1.0/sd_x + std::pow(x - mu_x, 2)/std::pow(sd_x, 3)) *
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
      stderror[i] = std::sqrt(avar(i,i)/(double)nobs);

  }


public:

  EmGaussGaussModel(
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
    integ_tol0 = integ_tol;
    integ_depth = integ_depth_;
    integ_method = integ_method_;

    // error catcher
    integ_tol_limit = 1e-14;
    accuracy_limit_error = false;

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
    sd_x = std::sqrt(s2 * 0.75);
    sigma = std::sqrt(s2 * 0.25);
    sd_w = std::sqrt(s2);

    // initialize avar and se
    avar = NumericMatrix(3, 3);
    rownames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    colnames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    stderror = NumericVector::create(
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
        Rcout << "iter " << i+1 << " : " <<
          " sigma = " << sigma << ", sd_x = " << sd_x <<
          " value = " << cur_value << "\n";
      }
      UpdateParameters();
      double increment = UpdateValueAndWeights();
      if (accuracy_limit_error) {
        // terminate the iteration if the accuracy limit has reached
        // did not converge
        warning("dit not converge: accuracy limit has reached");
        break;
      }

      if (std::fabs(increment) < tol*(std::fabs(cur_value - increment) + tol)) {
        convergence = 0;
        last_increment = increment;
        if (verbose) Rcout << "CONVERGED!\n";
        break;
      }
    }

    // use original integration accuracy for avar computation
    integ_tol = integ_tol0;
    UpdateAvarAndSe();
  }

  List CompileOutput()
  {
    List out;
    NumericVector estimate = NumericVector::create(
      Named("sigma") = sigma, Named("mu_x") = mu_x, Named("sd_x") = sd_x
    );

    out = List::create(
      Named("estimate") = estimate,
      Named("stderr") = stderror,
      Named("avar") = avar,
      Named("nobs") = nobs,
      Named("convergence") = convergence,
      Named("last_increment") = last_increment
    );
    return out;
  }

};


// [[Rcpp::export]]
List em_gauss_gauss_helper(
  std::vector<int> d_vec, std::vector<double> w_vec, double cutoff,
  double tol, int maxit,
  std::string integ_method, double integ_tol, int integ_depth,
  bool verbose)
{
  EmGaussGaussModel model(d_vec, w_vec, cutoff, tol, maxit,
                          integ_method, integ_tol, integ_depth);
  model.Estimate(verbose);
  return model.CompileOutput();
}


/*** R
library(rddsigma)
dat <- gen_data(500, 0.3, 1)
rddsigma:::em_gauss_gauss_helper(dat$d, dat$w, 1,
                               1e-6, 1000, "romberg", 1e-6, 100, TRUE)
rddsigma:::em_gauss_gauss(dat$d, dat$w, 1, reltol = 1e-6,
             integrate_options = list(rel.tol = 1e-6), quiet = FALSE)
*/
