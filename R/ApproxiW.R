#' ApproxiW
#' to be documented
#' @usage ApproxiW(WW, la, order)
#' @param WW  to be documented
#' @param la  to be documented
#' @param order  to be documented
#' @noRd
#' @return to be documented
ApproxiW <-
function (WW, la, order)
.Call("ApproxiW", WW, la, order, PACKAGE = "mgwrsar")
