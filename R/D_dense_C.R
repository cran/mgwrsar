#' D_dense_C
#' to be documented
#' @usage D_dense_C(xxi, yyi, xx, yy)
#' @param xxi  to be documented
#' @param yyi  to be documented
#' @param xx  to be documented
#' @param yy  to be documented
#' @keywords internal
#' @return to be documented
D_dense_C <-
function (xxi, yyi, xx, yy)
.Call("D_dense_C", xxi, yyi, xx, yy, PACKAGE = "mgwrsar")
