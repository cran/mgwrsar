#' C++ wrapper for MGWR mixed beta estimation via pivoted QR decomposition
#' @usage mgwr_beta_pivotal_qrp_mixed_cpp(XV, y, XC, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE, get_se=FALSE)
#' @param XV Varying coefficient design matrix.
#' @param y Response vector.
#' @param XC Fixed coefficient design matrix.
#' @param indexG Neighbor index matrix.
#' @param Wd Weight matrix.
#' @param TP Target points vector.
#' @param get_ts Logical, compute trace of hat matrix.
#' @param get_s Logical, compute hat matrix.
#' @param get_Rk Logical, compute Rk matrix.
#' @param get_se Logical, compute standard errors.
#' @keywords internal
#' @return A list with estimated coefficients and diagnostics.
#' @noRd
mgwr_beta_pivotal_qrp_mixed_cpp <- function(XV, y, XC, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE, get_se=FALSE) {
  .Call("_mgwrsar_gwr_beta_pivotal_qrp_cpp", XV, y, XC, indexG, Wd, TP, get_ts, get_s, get_Rk, get_se, PACKAGE = "mgwrsar")
}
