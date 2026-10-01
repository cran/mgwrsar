#' C++ wrapper for INST computation
#' @usage INST_C(A, B, C, D)
#' @param A,B,C,D Numeric matrices for the INST computation.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
INST_C <- function(A, B, C, D, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_INST_C", A, B, C, D, PACKAGE = "mgwrsar")
}
