#' Proj_C
#' to be documented
#' @usage Proj_C(HH, XX)
#' @param HH  to be documented
#' @param XX  to be documented
#' @keywords internal
#' @return to be documented
Proj_C <-
function (HH, XX)
.Call("Proj_C", HH, XX, PACKAGE = "mgwrsar")
