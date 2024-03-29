// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// omitNaCpp
NumericVector omitNaCpp(NumericVector x);
RcppExport SEXP _perturbLM_omitNaCpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(omitNaCpp(x));
    return rcpp_result_gen;
END_RCPP
}
// sortCpp
NumericVector sortCpp(NumericVector v);
RcppExport SEXP _perturbLM_sortCpp(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(sortCpp(v));
    return rcpp_result_gen;
END_RCPP
}
// calcPvalLessCpp
double calcPvalLessCpp(NumericVector v, double x);
RcppExport SEXP _perturbLM_calcPvalLessCpp(SEXP vSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPvalLessCpp(v, x));
    return rcpp_result_gen;
END_RCPP
}
// calcPvalGreaterCpp
double calcPvalGreaterCpp(NumericVector v, double x);
RcppExport SEXP _perturbLM_calcPvalGreaterCpp(SEXP vSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPvalGreaterCpp(v, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_perturbLM_omitNaCpp", (DL_FUNC) &_perturbLM_omitNaCpp, 1},
    {"_perturbLM_sortCpp", (DL_FUNC) &_perturbLM_sortCpp, 1},
    {"_perturbLM_calcPvalLessCpp", (DL_FUNC) &_perturbLM_calcPvalLessCpp, 2},
    {"_perturbLM_calcPvalGreaterCpp", (DL_FUNC) &_perturbLM_calcPvalGreaterCpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_perturbLM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
