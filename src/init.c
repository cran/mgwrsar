// ============================================================================
// init.c — Initialisation
// ============================================================================

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // pour NULL
#include <R_ext/Rdynload.h>

// ---------------------------------------------------------------------------
// 1. RcppExports_arma.cpp
// ---------------------------------------------------------------------------
extern SEXP _mgwrsar_gwr_beta_univar_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_gwr_beta_pivotal_qrp_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_mgwr_beta_pivotal_qrp_mixed_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// ---------------------------------------------------------------------------
// 2.  RcppExports_eigen.cpp
// ---------------------------------------------------------------------------
extern SEXP _mgwrsar_Proj_C(SEXP, SEXP);
extern SEXP _mgwrsar_Sl_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_INST_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_PhWY_C(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_QRcpp2_C(SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_ApproxiW(SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_mod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mgwrsar_knn_stable_sort(SEXP, SEXP);
extern SEXP _mgwrsar_compute_DS_DT_cpp(SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
  // Armadillo
  {"_mgwrsar_gwr_beta_univar_cpp",              (DL_FUNC) &_mgwrsar_gwr_beta_univar_cpp,               8},
  {"_mgwrsar_gwr_beta_pivotal_qrp_cpp",         (DL_FUNC) &_mgwrsar_gwr_beta_pivotal_qrp_cpp,          10},
  {"_mgwrsar_mgwr_beta_pivotal_qrp_mixed_cpp",  (DL_FUNC) &_mgwrsar_mgwr_beta_pivotal_qrp_mixed_cpp,   10},

  // Eigen
  {"_mgwrsar_Proj_C",                           (DL_FUNC) &_mgwrsar_Proj_C,                            2},
  {"_mgwrsar_Sl_C",                             (DL_FUNC) &_mgwrsar_Sl_C,                              4},
  {"_mgwrsar_INST_C",                           (DL_FUNC) &_mgwrsar_INST_C,                            4},
  {"_mgwrsar_PhWY_C",                           (DL_FUNC) &_mgwrsar_PhWY_C,                            4},
  {"_mgwrsar_QRcpp2_C",                         (DL_FUNC) &_mgwrsar_QRcpp2_C,                          3},
  {"_mgwrsar_ApproxiW",                         (DL_FUNC) &_mgwrsar_ApproxiW,                          3},
  {"_mgwrsar_mod",                              (DL_FUNC) &_mgwrsar_mod,                               10},
  {"_mgwrsar_knn_stable_sort",                  (DL_FUNC) &_mgwrsar_knn_stable_sort,                   2},
  {"_mgwrsar_compute_DS_DT_cpp",                (DL_FUNC) &_mgwrsar_compute_DS_DT_cpp,                 5},

  {NULL, NULL, 0}
};
void R_init_mgwrsar(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
