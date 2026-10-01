#' C++ wrapper for univariate GWR beta estimation
#' @usage gwr_beta_univar_cpp(y, x, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_se=FALSE, nthreads=1)
#' @param y Response vector.
#' @param x Predictor vector.
#' @param XV Varying coefficient design matrix.
#' @param indexG Neighbor index matrix.
#' @param Wd Weight matrix.
#' @param TP Target points vector.
#' @param get_ts Logical, compute trace of hat matrix.
#' @param get_s Logical, compute hat matrix.
#' @param get_se Logical, compute standard errors.
#' @param nthreads Integer; requested native thread count. Use 1 for sequential fallback.
#' @keywords internal
#' @return A list with estimated coefficients and diagnostics.
#' @noRd
gwr_beta_univar_cpp <- function(y, x, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE,get_se=FALSE, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_gwr_beta_univar_cpp", y, x, XV, indexG, Wd, TP, get_ts, get_s,get_se, PACKAGE = "mgwrsar")
}
