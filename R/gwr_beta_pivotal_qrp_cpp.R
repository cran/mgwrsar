#' to be documented
#' @usage gwr_beta_pivotal_qrp_cpp(X, y, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE, get_se=FALSE)
#' @keywords internal
#' @return to be documented
#' @noRd
gwr_beta_pivotal_qrp_cpp <- function(X, y, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE,get_se=FALSE) {
  .Call("_mgwrsar_gwr_beta_pivotal_qrp_cpp", X, y, XV, indexG, Wd, TP, get_ts, get_s, get_Rk, get_se, PACKAGE = "mgwrsar")
}
