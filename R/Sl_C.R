#' to be documented
#' @usage Sl_C(A, B, C, D)
#' @keywords internal
#' @return to be documented
#' @noRd
Sl_C <- function(A, B, C, D) {
  .Call("_mgwrsar_Sl_C", A, B, C, D, PACKAGE = "mgwrsar")
}
