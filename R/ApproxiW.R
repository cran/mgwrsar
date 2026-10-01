#' to be documented
#' @usage ApproxiW(A, B, C)
#' @keywords internal
#' @return to be documented
#' @noRd
ApproxiW <- function(W, TP, n, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_ApproxiW", W, TP, n, PACKAGE = "mgwrsar")
}
