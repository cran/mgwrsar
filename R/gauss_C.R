#' gauss_C
#' to be documented
#' @usage gauss_C(dd, hh, Minv)
#' @param dd  to be documented
#' @param hh  to be documented
#' @param Minv  to be documented
#' @keywords internal
#' @return to be documented
gauss_C <-
function (dd, hh, Minv)
.Call("gauss_C", dd, hh, Minv, PACKAGE = "mgwrsar")
