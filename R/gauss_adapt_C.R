#' Adaptive gaussian kernel, RcppEigen version
#' @usage gauss_adapt_C(d, h)
#' @param d  a vector of distance
#' @param h  bandwidth size expressed in number of neighbors
#' @return a vector of weight
#' @examples
#' \donttest{
#' w=gauss_adapt_C(-30:30,20)
#' plot(-30:30,w,type='l')
#' abline(v=-10)
#' abline(v=10)
#' }
gauss_adapt_C <-
function (d, h)
.Call("gauss_adapt_C", d, h, PACKAGE = "mgwrsar")
