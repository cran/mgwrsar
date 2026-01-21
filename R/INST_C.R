#' to be documented
#' @usage INST_C(A, B, C, D)
#' @keywords internal
#' @return to be documented
#' @noRd
INST_C <- function(A, B, C, D) {
  .Call("_mgwrsar_INST_C", A, B, C, D, PACKAGE = "mgwrsar")
}
