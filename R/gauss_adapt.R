#' Adaptive gaussian kernel
#' @usage gauss_adapt(d, h)
#' @param d  a vector of distance
#' @param h  bandwidth size expressed in number of neighbors
#' @return a vector of weight
#' @examples
#' \donttest{
#' w=gauss_adapt(-30:30,20)
#' plot(-30:30,w,type='l')
#' abline(v=-10)
#' abline(v=10)
#' }
gauss_adapt<-function(d,h){exp(-0.5 * (d/h)^2)}
