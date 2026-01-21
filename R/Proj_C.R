#' to be documented
#' @usage Proj_C(A, B)
#' @keywords internal
#' @return to be documented
#' @noRd
Proj_C <- function(A, B) {
  .Call("_mgwrsar_Proj_C", A, B, PACKAGE = "mgwrsar")
}
