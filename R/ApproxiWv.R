#' ApproxiWv
#' to be documented
#' @usage ApproxiWv(WW, la, order)
#' @param WW  to be documented
#' @param la  to be documented
#' @param order  to be documented
#' @noRd
#' @return to be documented
ApproxiWv <-
  function (WW, la, order)
    .Call("ApproxiWv", WW, la, order, PACKAGE = "mgwrsar")
