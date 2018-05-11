#' Sl_C
#' to be documented
#' @usage Sl_C(llambda, WW, iinv, aapprox)
#' @param llambda  to be documented
#' @param WW  to be documented
#' @param iinv  to be documented
#' @param aapprox  to be documented
#' @keywords internal
#' @return to be documented
Sl_C <-
function (llambda, WW, iinv, aapprox)
.Call("Sl_C", llambda, WW, iinv, aapprox, PACKAGE = "mgwrsar")
