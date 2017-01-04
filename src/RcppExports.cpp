// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// em_gauss_lap_helper
List em_gauss_lap_helper(std::vector<int> d_vec, std::vector<double> w_vec, double cutoff, double tol, int maxit, std::string integ_method, double integ_tol, int integ_depth, bool verbose);
RcppExport SEXP rddsigma_em_gauss_lap_helper(SEXP d_vecSEXP, SEXP w_vecSEXP, SEXP cutoffSEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP integ_methodSEXP, SEXP integ_tolSEXP, SEXP integ_depthSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type d_vec(d_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type w_vec(w_vecSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< std::string >::type integ_method(integ_methodSEXP);
    Rcpp::traits::input_parameter< double >::type integ_tol(integ_tolSEXP);
    Rcpp::traits::input_parameter< int >::type integ_depth(integ_depthSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(em_gauss_lap_helper(d_vec, w_vec, cutoff, tol, maxit, integ_method, integ_tol, integ_depth, verbose));
    return rcpp_result_gen;
END_RCPP
}
// integrate_test
void integrate_test();
RcppExport SEXP rddsigma_integrate_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    integrate_test();
    return R_NilValue;
END_RCPP
}
// integrate_test2
void integrate_test2(std::string method, double tol);
RcppExport SEXP rddsigma_integrate_test2(SEXP methodSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    integrate_test2(method, tol);
    return R_NilValue;
END_RCPP
}
// integrate_test3
void integrate_test3();
RcppExport SEXP rddsigma_integrate_test3() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    integrate_test3();
    return R_NilValue;
END_RCPP
}
// integrate_test4
void integrate_test4();
RcppExport SEXP rddsigma_integrate_test4() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    integrate_test4();
    return R_NilValue;
END_RCPP
}
// tsgauss_lfunc_helper
double tsgauss_lfunc_helper(IntegerVector d, NumericVector w, double cutoff, NumericVector mu_u, double sd_u, double sigma);
RcppExport SEXP rddsigma_tsgauss_lfunc_helper(SEXP dSEXP, SEXP wSEXP, SEXP cutoffSEXP, SEXP mu_uSEXP, SEXP sd_uSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_u(mu_uSEXP);
    Rcpp::traits::input_parameter< double >::type sd_u(sd_uSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(tsgauss_lfunc_helper(d, w, cutoff, mu_u, sd_u, sigma));
    return rcpp_result_gen;
END_RCPP
}
// tsgauss_lfunc_each_helper
NumericVector tsgauss_lfunc_each_helper(IntegerVector d, NumericVector w, double cutoff, NumericVector mu_u, double sd_u, double sigma);
RcppExport SEXP rddsigma_tsgauss_lfunc_each_helper(SEXP dSEXP, SEXP wSEXP, SEXP cutoffSEXP, SEXP mu_uSEXP, SEXP sd_uSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_u(mu_uSEXP);
    Rcpp::traits::input_parameter< double >::type sd_u(sd_uSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(tsgauss_lfunc_each_helper(d, w, cutoff, mu_u, sd_u, sigma));
    return rcpp_result_gen;
END_RCPP
}
