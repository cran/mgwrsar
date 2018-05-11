#' GDTC
#' to be documented
#' @usage GDTC(j, Y, X, S, H, kernels, minv = 0, maxknn = 500, NmaxDist = 6000,isgcv, TIME, decay)
#' @param j   to be documented
#' @param Y  to be documented
#' @param X  to be documented
#' @param S  to be documented
#' @param H  to be documented
#' @param kernels  to be documented
#' @param minv  to be documented
#' @param maxknn  to be documented
#' @param NmaxDist  to be documented
#' @param isgcv  to be documented
#' @param TIME  to be documented
#' @param decay  to be documented
#' @keywords internal
#'
#' @return to be documented
GDTC <-
function (j, Y, X, S, H, kernels, minv = 0, maxknn = 500, NmaxDist = 6000,
    isgcv, TIME, decay)
{
    Wz <- ifelse(S[, 4] != S[j + 1, 4], H[2 + S[j + 1, 4]], 1)
    kernel_C(as.matrix(S[, 1:2]), j, H[1], kernels[1], minv,
        FALSE, decay, isgcv, TRUE) * kernel_C(as.matrix(S[, 3]),
        j, H[2], kernels[2], minv, FALSE, decay, isgcv, TRUE) *
        Wz
}
