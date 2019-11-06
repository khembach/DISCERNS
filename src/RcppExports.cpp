// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_unique_ind
LogicalVector get_unique_ind(DataFrame& df);
RcppExport SEXP _exondiscovery_get_unique_ind(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(get_unique_ind(df));
    return rcpp_result_gen;
END_RCPP
}
// filter_exon_length
LogicalVector filter_exon_length(DataFrame& df, int max_length);
RcppExport SEXP _exondiscovery_filter_exon_length(SEXP dfSEXP, SEXP max_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type max_length(max_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(filter_exon_length(df, max_length));
    return rcpp_result_gen;
END_RCPP
}
// pred_exon_coord
DataFrame pred_exon_coord(DataFrame& df);
RcppExport SEXP _exondiscovery_pred_exon_coord(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_exon_coord(df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_exondiscovery_get_unique_ind", (DL_FUNC) &_exondiscovery_get_unique_ind, 1},
    {"_exondiscovery_filter_exon_length", (DL_FUNC) &_exondiscovery_filter_exon_length, 2},
    {"_exondiscovery_pred_exon_coord", (DL_FUNC) &_exondiscovery_pred_exon_coord, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_exondiscovery(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}