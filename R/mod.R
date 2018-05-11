#' mod
#' to be documented
#' @usage mod(YY, XX, WW, XZZ, YZZ, Wi, LocalInst, ismethodB2SLS, ismethodMGWRSAR_1_kc_0, SE_)
#' @param YY  to be documented
#' @param XX  to be documented
#' @param WW  to be documented
#' @param XZZ  to be documented
#' @param YZZ  to be documented
#' @param Wi  to be documented
#' @param LocalInst  to be documented
#' @param ismethodB2SLS  to be documented
#' @param ismethodMGWRSAR_1_kc_0  to be documented
#' @param SE_  to be documented
#' @keywords internal
#' @return to be documented
mod <-
function (YY, XX, WW, XZZ, YZZ, Wi, LocalInst, ismethodB2SLS,
    ismethodMGWRSAR_1_kc_0, SE_)
.Call("mod", YY, XX, WW, XZZ, YZZ, Wi, LocalInst, ismethodB2SLS,
    ismethodMGWRSAR_1_kc_0, SE_, PACKAGE = "mgwrsar")
