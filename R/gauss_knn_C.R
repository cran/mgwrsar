#' gauss_knn_C
#' to be documented
#' @usage gauss_knn_C(dd, hh)
#' @param dd  to be documented
#' @param hh  to be documented
#' @keywords internal
#' @return to be documented
gauss_knn_C <-
function (dd, hh)
.Call("gauss_knn_C", dd, hh, PACKAGE = "mgwrsar")
