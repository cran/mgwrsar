#' C++ wrapper for QR decomposition
#' @usage QRcpp2_C(A, B, C)
#' @param A,B,C Numeric matrices for the QR decomposition.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
QRcpp2_C <- function(A, B, C) {
  .Call("_mgwrsar_QRcpp2_C", A, B, C, PACKAGE = "mgwrsar")
}
