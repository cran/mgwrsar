#' PhWY_C
#' to be documented
#' @usage PhWY_C(YY, XX, WW, Wi)
#' @param YY  to be documented
#' @param XX  to be documented
#' @param WW  to be documented
#' @param Wi  to be documented
#' @keywords internal
#' @return to be documented
PhWY_C <-
function (YY, XX, WW, Wi)
.Call("PhWY_C", YY, XX, WW, Wi, PACKAGE = "mgwrsar")
