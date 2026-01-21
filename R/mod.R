#' to be documented
#' @usage mod(Y, X, W, X_s, Y_s, S, type, b2sls, get_ts, get_s)
#' @keywords internal
#' @return to be documented
#' @noRd
mod <- function(Y, X, W, X_s, Y_s, S, type, b2sls, get_ts, get_s) {
  .Call("_mgwrsar_mod", Y, X, W, X_s, Y_s, S, type, b2sls, get_ts, get_s, PACKAGE = "mgwrsar")
}
