#' C++ wrapper for PhWY computation
#' @usage PhWY_C(A, B, C, D)
#' @param A,B,C,D Numeric matrices for the PhWY computation.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
PhWY_C <- function(A, B, C, D) {
  .Call("_mgwrsar_PhWY_C", A, B, C, D, PACKAGE = "mgwrsar")
}
