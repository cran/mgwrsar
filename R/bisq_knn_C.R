#' Adaptive bisquare kernel
#' @usage bisq_knn_C(d,h)
#' @param d  a vector of distance
#' @param h  bandwidth size expressed in number of neighbors
#' @return a vector of weight
#' @examples
#' \donttest{
#' w=bisq_knn_C(-30:30,20)
#' plot(-30:30,w,type='l')
#' abline(v=-10)
#' abline(v=10)
#' }
bisq_knn_C <-
function (d, h)
.Call("bisq_knn_C", d, h, PACKAGE = "mgwrsar")
