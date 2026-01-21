#' to be documented
#' @usage mgwr_beta_pivotal_qrp_mixed_cpp(XV, y, XC, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE)
#' @keywords internal, get_se=FALSE
#' @return to be documented
#' @noRd
mgwr_beta_pivotal_qrp_mixed_cpp <- function(XV, y, XC, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE, get_Rk=FALSE, get_se=FALSE) {
  .Call("_mgwrsar_gwr_beta_pivotal_qrp_cpp", XV, y, XC, indexG, Wd, TP, get_ts, get_s, get_Rk, get_se, PACKAGE = "mgwrsar")
}
