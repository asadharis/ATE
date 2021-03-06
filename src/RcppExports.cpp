// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CRFamily
arma::vec CRFamily(arma::vec x, double theta);
RcppExport SEXP ATE_CRFamily(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(CRFamily(x, theta));
    return __result;
END_RCPP
}
// CRFamilyDerivative
arma::vec CRFamilyDerivative(arma::vec x, double theta);
RcppExport SEXP ATE_CRFamilyDerivative(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(CRFamilyDerivative(x, theta));
    return __result;
END_RCPP
}
// CRFamilySecondDerivative
arma::vec CRFamilySecondDerivative(arma::vec x, double theta);
RcppExport SEXP ATE_CRFamilySecondDerivative(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(CRFamilySecondDerivative(x, theta));
    return __result;
END_RCPP
}
// ObjectiveFunction
double ObjectiveFunction(NumericVector lam, NumericMatrix u, NumericVector ubar, arma::vec treat, double theta);
RcppExport SEXP ATE_ObjectiveFunction(SEXP lamSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP treatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(ObjectiveFunction(lam, u, ubar, treat, theta));
    return __result;
END_RCPP
}
// ObjectiveFirstDerivative
arma::rowvec ObjectiveFirstDerivative(NumericVector lam, NumericMatrix u, NumericVector ubar, arma::vec treat, double theta);
RcppExport SEXP ATE_ObjectiveFirstDerivative(SEXP lamSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP treatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(ObjectiveFirstDerivative(lam, u, ubar, treat, theta));
    return __result;
END_RCPP
}
// UpdateInverseHessian
arma::mat UpdateInverseHessian(arma::mat old_inv_hessian, arma::vec diff_in_est, arma::vec diff_in_derv);
RcppExport SEXP ATE_UpdateInverseHessian(SEXP old_inv_hessianSEXP, SEXP diff_in_estSEXP, SEXP diff_in_dervSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type old_inv_hessian(old_inv_hessianSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type diff_in_est(diff_in_estSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type diff_in_derv(diff_in_dervSEXP);
    __result = Rcpp::wrap(UpdateInverseHessian(old_inv_hessian, diff_in_est, diff_in_derv));
    return __result;
END_RCPP
}
// BacktrackLineSearch
double BacktrackLineSearch(double kAlpha, double kBeta, NumericVector current_est, NumericVector current_direction, NumericVector current_derv, NumericMatrix u, NumericVector ubar, arma::vec treat, double theta);
RcppExport SEXP ATE_BacktrackLineSearch(SEXP kAlphaSEXP, SEXP kBetaSEXP, SEXP current_estSEXP, SEXP current_directionSEXP, SEXP current_dervSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP treatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type kAlpha(kAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type kBeta(kBetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type current_est(current_estSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type current_direction(current_directionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type current_derv(current_dervSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(BacktrackLineSearch(kAlpha, kBeta, current_est, current_direction, current_derv, u, ubar, treat, theta));
    return __result;
END_RCPP
}
// BFGSAlgorithm
List BFGSAlgorithm(NumericVector initial, NumericMatrix u, NumericVector ubar, arma::vec treat, double theta, double kAlpha, double kBeta, int max_iter, double tol);
RcppExport SEXP ATE_BFGSAlgorithm(SEXP initialSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP treatSEXP, SEXP thetaSEXP, SEXP kAlphaSEXP, SEXP kBetaSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type treat(treatSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type kAlpha(kAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type kBeta(kBetaSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(BFGSAlgorithm(initial, u, ubar, treat, theta, kAlpha, kBeta, max_iter, tol));
    return __result;
END_RCPP
}
