#' C++ wrapper for projection computation
#' @usage Proj_C(A, B)
#' @param A,B Numeric matrices for the projection computation.
#' @keywords internal
#' @return A numeric matrix.
#' @noRd
Proj_C <- function(A, B) {
  .Call("_mgwrsar_Proj_C", A, B, PACKAGE = "mgwrsar")
}
