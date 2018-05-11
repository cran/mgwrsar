#' kernel_C
#' to be documented
#' @usage kernel_C(XX, J, hh, Mykernel, Minv, TIME, Decay, DDiagNull, normWW)
#' @param XX  to be documented
#' @param J  to be documented
#' @param hh  to be documented
#' @param Mykernel  to be documented
#' @param Minv  to be documented
#' @param TIME  to be documented
#' @param Decay  to be documented
#' @param DDiagNull  to be documented
#' @param normWW  to be documented
#' @keywords internal
#' @return to be documented
kernel_C <-
function (XX, J, hh, Mykernel, Minv, TIME, Decay, DDiagNull,
    normWW)
.Call("kernel_C", XX, J, hh, Mykernel, Minv, TIME, Decay, DDiagNull,
    normWW, PACKAGE = "mgwrsar")
