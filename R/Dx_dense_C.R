#' Dx_dense_C
#' to be documented
#' @usage Dx_dense_C(xxi, xx, TIME)
#' @param xxi  to be documented
#' @param xx  to be documented
#' @param TIME  to be documented
#' @keywords internal
#' @return to be documented
Dx_dense_C <-
function (xxi, xx, TIME)
.Call("Dx_dense_C", xxi, xx, TIME, PACKAGE = "mgwrsar")
