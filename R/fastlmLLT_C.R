#' fastlmLLT_C
#' to be documented
#' @usage fastlmLLT_C(XX, YY, SE_)
#' @param XX  to be documented
#' @param YY  to be documented
#' @param SE_  to be documented
#' @keywords internal
#'
#' @return to be documented
fastlmLLT_C <-
function (XX, YY, SE_)
.Call("fastlmLLT_C", XX, YY, SE_, PACKAGE = "mgwrsar")
