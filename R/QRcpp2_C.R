#' C++ wrapper for QR decomposition
#' @usage QRcpp2_C(A, B, C)
#' @param A,B,C Numeric matrices for the QR decomposition.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
QRcpp2_C <- function(A, B, C, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_QRcpp2_C", A, B, C, PACKAGE = "mgwrsar")
}
