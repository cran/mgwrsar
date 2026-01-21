#' to be documented
#' @usage ApproxiW(A, B, C)
#' @keywords internal
#' @return to be documented
#' @noRd
ApproxiW <- function(W, TP, n) {
  .Call("_mgwrsar_ApproxiW", W, TP, n, PACKAGE = "mgwrsar")
}
