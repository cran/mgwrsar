#' QRcpp2_C
#' to be documented
#' @usage QRcpp2_C(AA, bb, cc)
#' @param AA  to be documented
#' @param bb  to be documented
#' @param cc  to be documented
#' @keywords internal
#' @return to be documented
QRcpp2_C <-
function (AA, bb, cc)
.Call("QRcpp2_C", AA, bb, cc, PACKAGE = "mgwrsar")
