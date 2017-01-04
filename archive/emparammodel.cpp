#include <functional>
#include <cmath>
#include <Rcpp.h>
#include <string>
#include <vector>
#include "integral.h"
#include "emparammodel.h"
using namespace Rcpp;


// EmParamModel ------
void EmParamModel::Initialize(
    const std::vector<int> &d, const std::vector<double> &w, double c)
{
  d_vec = d;
  w_vec = w;
  cutoff = c;
  nobs = d_vec.size();

  if (w_vec.size() != d_vec.size()) stop("w and d must have the same length");
  if (nobs < 2) stop("too few data point");

  weights.resize(nobs, 1.0);
  cur_value = - INFINITY;

  convergence = -1;
}


void EmParamModel::Summary()
{
  Rcout << "Parameters:\n";
  for (int i = 0; i < theta_u.size(); i++)
    Rcout << " " << theta_u_names[i] << " = " << theta_u[i] << "\n";
  for (int i = 0; i < theta_x.size(); i++)
    Rcout << " " << theta_x_names[i] << " = " << theta_x[i] << "\n";
  Rcout << "Value = " << cur_value << "\n";
}


double EmParamModel::UpdateValueAndWeights(
    std::string integ_method, double integ_tol, int integ_depth)
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


void EmParamModel::UpdateParameters(
    std::string integ_method, double integ_tol, int integ_depth)
{
  std::vector<double> new_theta_u = GetNewThetaU(
    integ_method, integ_tol, integ_depth);
  std::vector<double> new_theta_x = GetNewThetaX(
    integ_method, integ_tol, integ_depth);

  theta_u = new_theta_u;
  theta_x = new_theta_x;
}


void EmParamModel::Estimate(
    double tol, int maxit, bool verbose,
    std::string integ_method, double integ_tol, int integ_depth)
{
  // update values and weights first using the initial parameters
  UpdateValueAndWeights(integ_method, integ_tol, integ_depth);

  // EM update
  convergence = 1;
  for (int i = 0; i < maxit; i++)
  {
    UpdateParameters(integ_method, integ_tol, integ_depth);
    double incr = UpdateValueAndWeights(integ_method, integ_tol, integ_depth);
    if (verbose) {
      Rcout << "iter " << i << "\n";
      Summary();
    }
    if (incr < tol*(cur_value - incr + tol)) {
      if (verbose) Rcout << "CONVERGED!\n";
      convergence = 0;
      break;
    }
  }
  Summary();

  // update asymptotic variance
  UpdateAvarAndSe();

}

// ------ EmParamModel




// EmParamGaussX -------
void EmParamGaussX::Initialize()
{
  //EmParamModel::Initialize(d, w, c);

  // initialize theta x
  theta_x.resize(2);
  theta_x_names.resize(2);
  theta_x_names[0] = "mu_x";
  theta_x_names[1] = "sd_x";
  double tmp1;
  double tmp2;
  for (int i = 0; i < nobs; i++)
  {
    tmp1 += w_vec[i];
    tmp2 += w_vec[i]*w_vec[i];
  }
  theta_x[0] = tmp1 / (double)nobs;
  theta_x[1] = sqrt(tmp2 / (double)nobs)*0.75;
}

double EmParamGaussX::fx(double x)
{
  return R::dnorm4(x, theta_x[0], theta_x[1], 0);
}

std::vector<double> EmParamGaussX::GetNewThetaX(
    std::string integ_method, double integ_tol, int integ_depth)
{
  double new_sdx = 0;
  for (int i = 0; i < nobs; i++)
  {
    std::function<double(double)> func = [this,i] (double x) -> double {
      return fx(x) * fu(x, i) / weights[i] * pow(x - theta_x[0], 2); };
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
  new_sdx = sqrt(new_sdx / (double)nobs);

  std::vector<double> out(2);
  out[0] = theta_x[0];
  out[1] = new_sdx;
  return out;
}
// ------ EmParamGaussX



// EmParamLapU -------
void EmParamLapU::Initialize()
{
  theta_u.resize(1);
  theta_u_names.resize(1);
  theta_u_names[0] = "sigma";

  double tmp;
  for (int i = 0; i < nobs; i++)
  {
    tmp += w_vec[i]*w_vec[i];
  }
  theta_u[0] = sqrt(tmp / (double)nobs)*0.25;
}


std::vector<double> EmParamLapU::GetNewThetaU(
    std::string integ_method, double integ_tol, int integ_depth)
{
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

  std::vector<double> out(1, new_sigma);
  return out;
}

// ------- EmParamLapU




// EmParamGaussLap -------
//
// EmParamGaussLap::EmParamGaussLap(
//     const std::vector<int> &d, const std::vector<double> &w, double c)
// {
//   // // store data
//   // EmParamModel::Initialize(d, w, c);
//   //
//   // // initialize parameters
//   // EmParamGaussX::Initialize();
//   // EmParamLapU::Initialize();
//
//   // // initialize avar and se
//   // Avar = NumericMatrix(3, 3);
//   // rownames(Avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
//   // colnames(Avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
//   //
//   // stderr = NumericVector::create(
//   //   Named("sigma") = 0,
//   //   Named("mu_x") = 0,
//   //   Named("sd_x") = 0);
// }

// void EmParamGaussLap::UpdateAvarAndSe()
// {
//   // todo: deploy
// }
// ------- EmParamGaussLap




// [[Rcpp::export]]
void em_model_test(
    std::vector<int> d_vec, std::vector<double> w_vec, double cutoff)
{
  double tol = 1e-5;
  int maxit = 1000;
  bool verbose = true;
  std::string integ_method = "romberg";
  double integ_tol = 1e-5;
  int integ_depth = 100;

  EmParamGaussLap model(d_vec, w_vec, cutoff);
  // model.Summary();
  // //model.Estimate(tol, maxit, verbose, integ_method, integ_tol, integ_depth);
  // model.Summary();
}

/*** R
library(rddsigma)
dat <- gen_data(500, 0.3, 1)
rddsigma:::em_model_test(dat$d, dat$w, 1)
*/