#' EDK
#' to be documented
#' @usage EDK(d, h, k = 2)
#' @param d  to be documented
#' @param h  to be documented
#' @param k  to be documented
#' @keywords internal
#' @return to be documented
EDK <-
function (d, h, k = 2)
{
    exp(-0.5 * (d/h)^k)
}
