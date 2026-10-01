#' C++ wrapper for Sl computation
#' @usage Sl_C(A, B, C, D)
#' @param A,B,C,D Numeric matrices for the Sl computation.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
Sl_C <- function(A, B, C, D, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_Sl_C", A, B, C, D, PACKAGE = "mgwrsar")
}
