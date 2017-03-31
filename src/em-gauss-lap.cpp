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
  { return exp(-std::sqrt(2.0) / sigma * std::fabs(w_vec[i] - x)) / sigma / std::sqrt(2.0); }

  double new_fx(double x)
  { return R::dnorm4(x, mu_x, new_sdx, 0); }
  double new_fu(double x, int i)
  { return exp(-std::sqrt(2.0) / new_sigma * std::fabs(w_vec[i] - x)) / new_sigma / std::sqrt(2.0); }

  // temporary storage of updated params, value, weights
  double new_sigma;
  double new_sdx;
  double new_value;
  std::vector<double> new_weights;

  // momentum
  double mom_sdx;
  double mom_sigma;



  double Update(bool verbose)
  {
    // Update (param, weights, value)
    // The process make sures that the value is improved
    // (up to small numerical error)
    // returns the increments in the value

    ComputeNewParameters();
    ComputeNewValueAndWeights();

    double increment;
    //integ_tol = integ_tol0;
    while (true)
    {
      increment = new_value - cur_value;
      // is increment positive (up to some small error)?
      // theoretically EM update must improve the value, but
      // numerical problem may occur
      if (increment > -tol*(std::fabs(cur_value) + tol)) break;


      // First, try increasing integration accuracy up to the limit
      if (integ_tol < integ_tol_limit) break;
      integ_tol *= 0.2;
      if (verbose) {
        Rcout.precision(5);
        Rcout << "negative increment: " << increment;
        Rcout << " integ_tol is now " << integ_tol << "\n";
        Rcout.flush();
      }
      // first, recompute value and weights with the new integ_tol level
      RecomputeValueAndWeights();
      ComputeNewValueAndWeights();
    }

    double mom_rate = 1;
    int count = 0;
    while (true)
    {
      increment = new_value - cur_value;
      if (increment > -tol*(std::fabs(cur_value) + tol)) break;

      // as the last resort, we change the parameters
      // around the current values to see if value gets improved
      // magnitude of change (mom_rate) starts from 1
      // and gradually becomes smaller
      // Eventually the rate is almost zero, i.e. no paramter change
      if (verbose) {
        Rcout.precision(5);
        Rcout << "negative increment: " << increment;
        Rcout << " seek the neighborhood with rate " << mom_rate << "\n";
        //Rcout << "sigma = " << new_sigma << ", value = " << new_value << "\n";
        Rcout.flush();
      }
      int sign1 = count % 2 == 0 ? 1 : -1;
      new_sigma = sigma + sign1*mom_rate * mom_sigma;
      int sign2 = count % 4 < 2 ? 1 : -1;
      new_sdx   = sd_x + sign2*mom_rate * mom_sdx;

      count++;
      if (count >= 4) {
        mom_rate *= 0.8;
        count = 0;
      }

      ComputeNewValueAndWeights();
    }


    // now the value has been improved.
    // so apply new (params, weights, value)
    mom_sigma = new_sigma - sigma;
    mom_sdx = new_sdx - sd_x;
    sigma = new_sigma;
    sd_x = new_sdx;
    cur_value = new_value;
    for (size_t i = 0; i < weights.size(); i++) weights[i] = new_weights[i];

    return increment;
  }



  void RecomputeValueAndWeights()
  {
    cur_value = 0;

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
      cur_value += log(weights[i]);
    }
    cur_value /= double (nobs);
  }


  void ComputeNewValueAndWeights()
  {
    new_value = 0;

    for (int i = 0; i < nobs; i++)
    {
      std::function<double(double)> func = [this, i] (double x) -> double {
        return new_fx(x) * new_fu(x, i); };

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
  }

  void ComputeNewParameters()
  {
    // update parameters for u
    new_sigma = 0;
    for (int i = 0; i < nobs; i++)
    {
      std::function<double(double)> func = [this,i] (double x) -> double {
        return fx(x) * fu(x, i) / weights[i] * std::fabs(w_vec[i] - x); };
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
    new_sigma *= (std::sqrt(2.0) / (double)nobs);

    // update paramters for x
    new_sdx = 0;
    for (int i = 0; i < nobs; i++)
    {
      std::function<double(double)> func = [this,i] (double x) -> double {
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
        return (-1.0/sigma + std::sqrt(2)*std::fabs(w_vec[i]-x)/std::pow(sigma, 2)) *
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
      stderror[i] = std::sqrt(avar(i,i)/nobs);

  }


public:

  EmGaussLapModel(
    const std::vector<int> &d_vec_, const std::vector<double> &w_vec_,
    double cutoff_,
    double init_sigma,
    double tol_, int maxit_,
    std::string integ_method_, double integ_tol_, int integ_depth_)
  {
    d_vec = d_vec_;
    w_vec = w_vec_;
    nobs = d_vec.size();
    cutoff = cutoff_;
    tol = tol_;
    maxit = maxit_,
    integ_tol = integ_tol_;
    integ_tol0 = integ_tol_;
    integ_depth = integ_depth_;
    integ_method = integ_method_;

    // error catcher
    integ_tol_limit = 1e-10;
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
    //sd_x = std::sqrt(s2 * 0.75);
    //sigma = std::sqrt(s2 * 0.25);
    sd_w = std::sqrt(s2);
    sigma = init_sigma;
    sd_x = std::sqrt(s2 - sigma*sigma);

    // initialize avar and se
    avar = NumericMatrix(3, 3);
    rownames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    colnames(avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
    stderror = NumericVector::create(
      Named("sigma") = 0, Named("mu_x") = 0, Named("sd_x") = 0);


    // initialize value and weights
    cur_value = -INFINITY;
    weights.resize(nobs, 1);

    new_sdx = sd_x;
    new_sigma = sigma;
    new_weights.resize(nobs, 1);
    new_value = -INFINITY;

    mom_sdx = sd_x*0.1;
    mom_sigma = sigma*0.1;

    convergence = -1;
  }

  void Estimate(bool verbose)
  {
    // initialize value and weights with the initial params
    RecomputeValueAndWeights();

    convergence = 1;
    for (int i = 0; i < maxit; i++)
    {
      if (verbose) {
        Rcout.precision(3);
        Rcout << "iter " << i+1 << " : " <<
          " sigma = " << sigma << ", sd_x = " << sd_x;
        Rcout.precision(5);
        Rcout << " value = " << cur_value << "\n";
      }
      double increment = Update(verbose);

      if (accuracy_limit_error) {
        // terminate the iteration if the accuracy limit has reached
        // did not converge
        warning("dit not converge: reached the accuracy limit");
        break;
      }

      if (std::fabs(increment) < tol*(std::fabs(cur_value - increment) + tol)) {
        convergence = 0;
        if (verbose) Rcout << "CONVERGED!\n";
        break;
      }
    }

    // use original accuracy level for computing avar
    integ_tol = integ_tol0;
    UpdateAvarAndSe();
  }

  List CompileOutput()
  {
    NumericVector estimate = NumericVector::create(
      Named("sigma") = sigma, Named("mu_x") = mu_x, Named("sd_x") = sd_x
    );

    List out = List::create(
      Named("estimate") = estimate,
      Named("stderr") = stderror,
      Named("avar") = avar,
      Named("nobs") = nobs,
      Named("convergence") = convergence,
      Named("value") = cur_value,
      Named("last_increment") = last_increment
    );
    return out;
  }

};


// [[Rcpp::export]]
List em_gauss_lap_helper(
    std::vector<int> d_vec, std::vector<double> w_vec, double cutoff,
    double init_sigma,
    double tol, int maxit,
    std::string integ_method, double integ_tol, int integ_depth,
    bool verbose)
{
  EmGaussLapModel model(d_vec, w_vec, cutoff,
                        init_sigma,
                        tol, maxit,
                        integ_method, integ_tol, integ_depth);
  model.Estimate(verbose);
  return model.CompileOutput();
}


/*** R
library(rddsigma)
dat <- gen_data(500, 0.2, 1)
rddsigma:::em_gauss_lap_helper(dat$d, dat$w, 1, 0.5,
                               1e-5, 1000, "romberg", 1e-6, 100, TRUE)
rddsigma:::em_gauss_lap_helper(dat$d, dat$w, 1, 0.1,
                               1e-5, 1000, "romberg", 1e-6, 100, TRUE)

dat <- gen_data(500, 0.2, 1, x_dist = "exp")
rddsigma:::em_gauss_lap_helper(dat$d, dat$w, 1, 0.5,
                               1e-5, 1000, "romberg", 1e-6, 100, TRUE)
rddsigma:::em_gauss_lap_helper(dat$d, dat$w, 1, 0.1,
                               1e-5, 1000, "romberg", 1e-6, 100, TRUE)

*/
