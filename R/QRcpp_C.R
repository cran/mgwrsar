#' QRcpp_C
#' to be documented
#' @usage QRcpp_C(AA, bb)
#' @param AA  to be documented
#' @param bb  to be documented
#' @keywords internal
#' @return to be documented
QRcpp_C <-
function (AA, bb)
.Call("QRcpp_C", AA, bb, PACKAGE = "mgwrsar")
