#' bin_C
#' to be documented
#' @usage bin_C(dd, hh, Minv)
#' @param dd  to be documented
#' @param hh  to be documented
#' @param Minv  to be documented
#' @keywords internal
#' @return to be documented
bin_C <-
function (dd, hh, Minv)
.Call("bin_C", dd, hh, Minv, PACKAGE = "mgwrsar")
