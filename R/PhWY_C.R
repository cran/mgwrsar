#' to be documented
#' @usage PhWY_C(A, B, C, D)
#' @keywords internal
#' @return to be documented
#' @noRd
PhWY_C <- function(A, B, C, D) {
  .Call("_mgwrsar_PhWY_C", A, B, C, D, PACKAGE = "mgwrsar")
}
