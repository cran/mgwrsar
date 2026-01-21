#' to be documented
#' @usage QRcpp2_C(A, B, C)
#' @keywords internal
#' @return to be documented
#' @noRd
QRcpp2_C <- function(A, B, C) {
  .Call("_mgwrsar_QRcpp2_C", A, B, C, PACKAGE = "mgwrsar")
}
