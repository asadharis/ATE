// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cpp_cr_rho
arma::vec cpp_cr_rho(arma::vec x, double theta);
RcppExport SEXP ATE_cpp_cr_rho(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_cr_rho(x, theta));
    return __result;
END_RCPP
}
// cpp_dcr_rho
arma::vec cpp_dcr_rho(arma::vec x, double theta);
RcppExport SEXP ATE_cpp_dcr_rho(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_dcr_rho(x, theta));
    return __result;
END_RCPP
}
// cpp_ddcr_rho
arma::vec cpp_ddcr_rho(arma::vec x, double theta);
RcppExport SEXP ATE_cpp_ddcr_rho(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_ddcr_rho(x, theta));
    return __result;
END_RCPP
}
// cpp_obj
double cpp_obj(NumericVector lam, NumericMatrix u, NumericVector ubar, arma::vec Ti, double theta);
RcppExport SEXP ATE_cpp_obj(SEXP lamSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP TiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ti(TiSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_obj(lam, u, ubar, Ti, theta));
    return __result;
END_RCPP
}
// cpp_derv_obj
arma::rowvec cpp_derv_obj(NumericVector lam, NumericMatrix u, NumericVector ubar, arma::vec Ti, double theta);
RcppExport SEXP ATE_cpp_derv_obj(SEXP lamSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP TiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ti(TiSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_derv_obj(lam, u, ubar, Ti, theta));
    return __result;
END_RCPP
}
// cpp_update_hessianInv
arma::mat cpp_update_hessianInv(arma::mat oldInv, arma::vec sk, arma::vec yk);
RcppExport SEXP ATE_cpp_update_hessianInv(SEXP oldInvSEXP, SEXP skSEXP, SEXP ykSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type oldInv(oldInvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sk(skSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yk(ykSEXP);
    __result = Rcpp::wrap(cpp_update_hessianInv(oldInv, sk, yk));
    return __result;
END_RCPP
}
// cpp_backtrack
double cpp_backtrack(double alpha, double beta, NumericVector x_val, NumericVector del_x, NumericVector nabla_f, NumericMatrix u, NumericVector ubar, arma::vec Ti, double theta);
RcppExport SEXP ATE_cpp_backtrack(SEXP alphaSEXP, SEXP betaSEXP, SEXP x_valSEXP, SEXP del_xSEXP, SEXP nabla_fSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP TiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_val(x_valSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del_x(del_xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nabla_f(nabla_fSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ti(TiSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    __result = Rcpp::wrap(cpp_backtrack(alpha, beta, x_val, del_x, nabla_f, u, ubar, Ti, theta));
    return __result;
END_RCPP
}
// cpp_quasi_newt
List cpp_quasi_newt(NumericVector ini, NumericMatrix u, NumericVector ubar, arma::vec Ti, double theta, double alpha, double beta, int max_iter, double tol);
RcppExport SEXP ATE_cpp_quasi_newt(SEXP iniSEXP, SEXP uSEXP, SEXP ubarSEXP, SEXP TiSEXP, SEXP thetaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type ini(iniSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ubar(ubarSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ti(TiSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(cpp_quasi_newt(ini, u, ubar, Ti, theta, alpha, beta, max_iter, tol));
    return __result;
END_RCPP
}
