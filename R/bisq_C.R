#' bisquare kernel, RcppEigen version
#' @usage bisq_C(d, h, Minv)
#' @param d  a vector of distance
#' @param h  a distance bandwidth
#' @param Minv  Minimal number of non null weight (neigbor)
#' @return a vector of weight
#' @examples
#' \donttest{
#' w=bisq_C(1:100,20,0)
#' plot(1:100,w,type='l')
#' }
bisq_C <-
function (d, h, Minv)
.Call("bisq_C", d, h, Minv, PACKAGE = "mgwrsar")
