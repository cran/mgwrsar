#' to be documented
#' @usage gwr_beta_univar_cpp(y, x, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE)
#' @keywords internal
#' @return to be documented
#' @noRd
gwr_beta_univar_cpp <- function(y, x, XV, indexG, Wd, TP, get_ts=FALSE, get_s=FALSE) {
  .Call("_mgwrsar_gwr_beta_univar_cpp", y, x, XV, indexG, Wd, TP, get_ts, get_s, PACKAGE = "mgwrsar")
}
