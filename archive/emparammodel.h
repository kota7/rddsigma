#include <functional>
#include <cmath>
#include <Rcpp.h>
#include <string>
#include <vector>
#include "integral.h"
using namespace Rcpp;



class EmParamModel
{

protected:

  // data
  std::vector<int> d_vec;
  std::vector<double> w_vec;
  int nobs;
  double cutoff;

  // model parameters
  std::vector<double> theta_x;
  std::vector<double> theta_u;
  std::vector<std::string> theta_x_names;
  std::vector<std::string> theta_u_names;

  // asymptoric variance for the parameters
  // and the standard errors
  NumericMatrix Avar;
  NumericVector stderr;

  // value and weights
  double cur_value;
  std::vector<double> weights;
  int convergence;

  // density functions as pure virtual,
  // to be defined in sub classes
  virtual double fx(double x) = 0;
  virtual double fu(double x, int i) = 0;


  // initializer, or data receiver
  void Initialize(
      const std::vector<int> &d, const std::vector<double> &w, double c);

  // update value and weights based on the current parameters
  // and returns the increment.
  // can define a generic function applicable for sub classes
  virtual double UpdateValueAndWeights(
      std::string integ_method, double integ_tol, int integ_depth);

  // parameter updaters.
  // 'get' functions compute the updated parameters, which depend on subclass
  // UpdateParamters routine can be generic
  virtual void UpdateParameters(
      std::string integ_method, double integ_tol, int integ_depth);
  virtual std::vector<double> GetNewThetaX(
      std::string integ_method, double integ_tol, int integ_depth) = 0;
  virtual std::vector<double> GetNewThetaU(
      std::string integ_method, double integ_tol, int integ_depth) = 0;

  // Asymptotic variance calculator based on current parameters
  virtual void UpdateAvarAndSe() = 0;


public:

  // EmParamModel(
  //   const std::vector<int> &d, const std::vector<double> &w, double c);
  virtual ~EmParamModel() {} ;

  // estimation routine, this can be generic to subclass
  virtual void Estimate(
      double tol, int maxit, bool verbose,
      std::string integ_method, double integ_tol, int integ_depth);

  // convenient function to display current status
  virtual void Summary();
};




// middle classes are defined for each x distribution and u distribution
// they are still abstruct.
class EmParamGaussX : public virtual EmParamModel
{
protected:
  double fx(double x);
  std::vector<double> GetNewThetaX(
      std::string integ_method, double integ_tol, int integ_depth);
  void Initialize();
};


class EmParamLapU : public virtual EmParamModel
{
protected:
  double fu(double x, int i);
  std::vector<double> GetNewThetaU(
      std::string integ_method, double integ_tol, int integ_depth);
  void Initialize();
};



// non-abstruct classes are defined for each combination of
// x and u distribution
class EmParamGaussLap : public EmParamGaussX, public EmParamLapU
{
// protected:
   void UpdateAvarAndSe() {}
public:
  EmParamGaussLap(
    const std::vector<int> &d, const std::vector<double> &w, double c)
  {
      EmParamModel::Initialize(d, w, c);

      // initialize parameters
      EmParamGaussX::Initialize();
      EmParamLapU::Initialize();

      // initialize avar and se
      Avar = NumericMatrix(3, 3);
      rownames(Avar) = CharacterVector::create("sigma", "mu_x", "sd_x");
      colnames(Avar) = CharacterVector::create("sigma", "mu_x", "sd_x");

      stderr = NumericVector::create(
        Named("sigma") = 0,
        Named("mu_x") = 0,
        Named("sd_x") = 0);
  }
  ~EmParamGaussLap() {}
};




