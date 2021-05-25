// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pgamma_ssd
double pgamma_ssd(double q, double shape, double scale);
RcppExport SEXP _ssdtools_pgamma_ssd(SEXP qSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(pgamma_ssd(q, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// qgamma_ssd
double qgamma_ssd(double p, double shape, double scale);
RcppExport SEXP _ssdtools_qgamma_ssd(SEXP pSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(qgamma_ssd(p, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// rgamma_ssd
NumericVector rgamma_ssd(int n, double shape, double scale);
RcppExport SEXP _ssdtools_rgamma_ssd(SEXP nSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rgamma_ssd(n, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// pgompertz_ssd
double pgompertz_ssd(double q, double location, double shape);
RcppExport SEXP _ssdtools_pgompertz_ssd(SEXP qSEXP, SEXP locationSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(pgompertz_ssd(q, location, shape));
    return rcpp_result_gen;
END_RCPP
}
// qgompertz_ssd
double qgompertz_ssd(double p, double location, double shape);
RcppExport SEXP _ssdtools_qgompertz_ssd(SEXP pSEXP, SEXP locationSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(qgompertz_ssd(p, location, shape));
    return rcpp_result_gen;
END_RCPP
}
// rgompertz_ssd
NumericVector rgompertz_ssd(int n, double location, double shape);
RcppExport SEXP _ssdtools_rgompertz_ssd(SEXP nSEXP, SEXP locationSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(rgompertz_ssd(n, location, shape));
    return rcpp_result_gen;
END_RCPP
}
// pgumbel_ssd
double pgumbel_ssd(double q, double location, double scale);
RcppExport SEXP _ssdtools_pgumbel_ssd(SEXP qSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(pgumbel_ssd(q, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// qgumbel_ssd
double qgumbel_ssd(double p, double location, double scale);
RcppExport SEXP _ssdtools_qgumbel_ssd(SEXP pSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(qgumbel_ssd(p, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// rgumbel_ssd
NumericVector rgumbel_ssd(int n, double location, double scale);
RcppExport SEXP _ssdtools_rgumbel_ssd(SEXP nSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rgumbel_ssd(n, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// plnorm_ssd
double plnorm_ssd(double q, double meanlog, double sdlog);
RcppExport SEXP _ssdtools_plnorm_ssd(SEXP qSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< double >::type sdlog(sdlogSEXP);
    rcpp_result_gen = Rcpp::wrap(plnorm_ssd(q, meanlog, sdlog));
    return rcpp_result_gen;
END_RCPP
}
// qlnorm_ssd
double qlnorm_ssd(double p, double meanlog, double sdlog);
RcppExport SEXP _ssdtools_qlnorm_ssd(SEXP pSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< double >::type sdlog(sdlogSEXP);
    rcpp_result_gen = Rcpp::wrap(qlnorm_ssd(p, meanlog, sdlog));
    return rcpp_result_gen;
END_RCPP
}
// rlnorm_ssd
NumericVector rlnorm_ssd(int n, double meanlog, double sdlog);
RcppExport SEXP _ssdtools_rlnorm_ssd(SEXP nSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< double >::type sdlog(sdlogSEXP);
    rcpp_result_gen = Rcpp::wrap(rlnorm_ssd(n, meanlog, sdlog));
    return rcpp_result_gen;
END_RCPP
}
// plogis_ssd
double plogis_ssd(double q, double location, double scale);
RcppExport SEXP _ssdtools_plogis_ssd(SEXP qSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(plogis_ssd(q, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// qlogis_ssd
double qlogis_ssd(double p, double location, double scale);
RcppExport SEXP _ssdtools_qlogis_ssd(SEXP pSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(qlogis_ssd(p, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// rlogis_ssd
NumericVector rlogis_ssd(int n, double location, double scale);
RcppExport SEXP _ssdtools_rlogis_ssd(SEXP nSEXP, SEXP locationSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rlogis_ssd(n, location, scale));
    return rcpp_result_gen;
END_RCPP
}
// pweibull_ssd
double pweibull_ssd(double q, double shape, double scale);
RcppExport SEXP _ssdtools_pweibull_ssd(SEXP qSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(pweibull_ssd(q, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// qweibull_ssd
double qweibull_ssd(double p, double shape, double scale);
RcppExport SEXP _ssdtools_qweibull_ssd(SEXP pSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(qweibull_ssd(p, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// rweibull_ssd
NumericVector rweibull_ssd(int n, double shape, double scale);
RcppExport SEXP _ssdtools_rweibull_ssd(SEXP nSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rweibull_ssd(n, shape, scale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ssdtools_pgamma_ssd", (DL_FUNC) &_ssdtools_pgamma_ssd, 3},
    {"_ssdtools_qgamma_ssd", (DL_FUNC) &_ssdtools_qgamma_ssd, 3},
    {"_ssdtools_rgamma_ssd", (DL_FUNC) &_ssdtools_rgamma_ssd, 3},
    {"_ssdtools_pgompertz_ssd", (DL_FUNC) &_ssdtools_pgompertz_ssd, 3},
    {"_ssdtools_qgompertz_ssd", (DL_FUNC) &_ssdtools_qgompertz_ssd, 3},
    {"_ssdtools_rgompertz_ssd", (DL_FUNC) &_ssdtools_rgompertz_ssd, 3},
    {"_ssdtools_pgumbel_ssd", (DL_FUNC) &_ssdtools_pgumbel_ssd, 3},
    {"_ssdtools_qgumbel_ssd", (DL_FUNC) &_ssdtools_qgumbel_ssd, 3},
    {"_ssdtools_rgumbel_ssd", (DL_FUNC) &_ssdtools_rgumbel_ssd, 3},
    {"_ssdtools_plnorm_ssd", (DL_FUNC) &_ssdtools_plnorm_ssd, 3},
    {"_ssdtools_qlnorm_ssd", (DL_FUNC) &_ssdtools_qlnorm_ssd, 3},
    {"_ssdtools_rlnorm_ssd", (DL_FUNC) &_ssdtools_rlnorm_ssd, 3},
    {"_ssdtools_plogis_ssd", (DL_FUNC) &_ssdtools_plogis_ssd, 3},
    {"_ssdtools_qlogis_ssd", (DL_FUNC) &_ssdtools_qlogis_ssd, 3},
    {"_ssdtools_rlogis_ssd", (DL_FUNC) &_ssdtools_rlogis_ssd, 3},
    {"_ssdtools_pweibull_ssd", (DL_FUNC) &_ssdtools_pweibull_ssd, 3},
    {"_ssdtools_qweibull_ssd", (DL_FUNC) &_ssdtools_qweibull_ssd, 3},
    {"_ssdtools_rweibull_ssd", (DL_FUNC) &_ssdtools_rweibull_ssd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ssdtools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
