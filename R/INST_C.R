#' INST_C
#' to be documented
#' @usage INST_C(XX, WW, withlambda, llambda)
#' @param XX  to be documented
#' @param WW  to be documented
#' @param withlambda  to be documented
#' @param llambda  to be documented
#' @keywords internal
#' @return to be documented
INST_C <-
function (XX, WW, withlambda, llambda)
.Call("INST_C", XX, WW, withlambda, llambda, PACKAGE = "mgwrsar")
