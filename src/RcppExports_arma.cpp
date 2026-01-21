// ============================================================================
// RcppExports_arma.cpp — Export RcppArmadillo routines for mgwrsar
// ============================================================================

#include "RcppIncludesArmadillo.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Declarations of exposed C++ functions (defined in gwr_core.cpp)
// ---------------------------------------------------------------------------

// C++ interfaces (with Rcpp types) defined in gwr_core.cpp
Rcpp::List gwr_beta_univar_cpp(const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector& x,
                               const Rcpp::NumericMatrix& XV,
                               const Rcpp::IntegerMatrix& indexG,
                               const Rcpp::NumericMatrix& Wd,
                               const Rcpp::IntegerVector& TP,
                               bool get_ts,
                               bool get_s);

Rcpp::List gwr_beta_pivotal_qrp_cpp(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& y,
                                    const Rcpp::NumericMatrix& XV,
                                    const Rcpp::IntegerMatrix& indexG,
                                    const Rcpp::NumericMatrix& Wd,
                                    const Rcpp::IntegerVector& TP,
                                    bool get_ts,
                                    bool get_s,
                                    bool get_Rk,
                                    bool get_se);

Rcpp::List mgwr_beta_pivotal_qrp_mixed_cpp(const Rcpp::NumericMatrix& XV,
                                           const Rcpp::NumericVector& y,
                                           const Rcpp::NumericMatrix& XC,
                                           const Rcpp::IntegerMatrix& indexG,
                                           const Rcpp::NumericMatrix& Wd,
                                           const Rcpp::IntegerVector& TP,
                                           bool get_ts,
                                           bool get_s,
                                           bool get_Rk,
                                           bool get_se);


// ---------------------------------------------------------------------------
// C wrappers (SEXP) for .Call() — names used by RcppExports.R
// ---------------------------------------------------------------------------

extern "C" {

  // gwr_beta_univar_cpp
  RcppExport SEXP _mgwrsar_gwr_beta_univar_cpp(SEXP ySEXP, SEXP xSEXP,
                                               SEXP XVSEXP, SEXP indexGSEXP,
                                               SEXP WdSEXP, SEXP TPSEXP,
                                               SEXP get_tsSEXP, SEXP get_sSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type XV(XVSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type indexG(indexGSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Wd(WdSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type TP(TPSEXP);
    Rcpp::traits::input_parameter< bool >::type get_ts(get_tsSEXP);
    Rcpp::traits::input_parameter< bool >::type get_s(get_sSEXP);
    rcpp_result_gen = Rcpp::wrap(
      gwr_beta_univar_cpp(y, x, XV, indexG, Wd, TP, get_ts, get_s)
    );
    return rcpp_result_gen;
    END_RCPP
  }

  // gwr_beta_pivotal_qrp_cpp
  RcppExport SEXP _mgwrsar_gwr_beta_pivotal_qrp_cpp(SEXP XSEXP, SEXP ySEXP,
                                                    SEXP XVSEXP, SEXP indexGSEXP,
                                                    SEXP WdSEXP, SEXP TPSEXP,
                                                    SEXP get_tsSEXP, SEXP get_sSEXP,
                                                    SEXP get_RkSEXP,SEXP get_seSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type XV(XVSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type indexG(indexGSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Wd(WdSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type TP(TPSEXP);
    Rcpp::traits::input_parameter< bool >::type get_ts(get_tsSEXP);
    Rcpp::traits::input_parameter< bool >::type get_s(get_sSEXP);
    Rcpp::traits::input_parameter< bool >::type get_Rk(get_RkSEXP);
    Rcpp::traits::input_parameter< bool >::type get_se(get_seSEXP);
    rcpp_result_gen = Rcpp::wrap(
      gwr_beta_pivotal_qrp_cpp(X, y, XV, indexG, Wd, TP, get_ts, get_s, get_Rk,get_se)
    );
    return rcpp_result_gen;
    END_RCPP
  }

  // mgwr_beta_pivotal_qrp_mixed_cpp
  RcppExport SEXP _mgwrsar_mgwr_beta_pivotal_qrp_mixed_cpp(SEXP XVSEXP, SEXP ySEXP,
                                                           SEXP XCSEXP, SEXP indexGSEXP,
                                                           SEXP WdSEXP, SEXP TPSEXP,
                                                           SEXP get_tsSEXP, SEXP get_sSEXP,
                                                           SEXP get_RkSEXP,SEXP get_seSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type XV(XVSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type XC(XCSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type indexG(indexGSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Wd(WdSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type TP(TPSEXP);
    Rcpp::traits::input_parameter< bool >::type get_ts(get_tsSEXP);
    Rcpp::traits::input_parameter< bool >::type get_s(get_sSEXP);
    Rcpp::traits::input_parameter< bool >::type get_Rk(get_RkSEXP);
    Rcpp::traits::input_parameter< bool >::type get_se(get_seSEXP);
    rcpp_result_gen = Rcpp::wrap(
      mgwr_beta_pivotal_qrp_mixed_cpp(XV, y, XC, indexG, Wd, TP, get_ts, get_s, get_Rk,get_se)
    );
    return rcpp_result_gen;
    END_RCPP
  }

} // extern "C"
