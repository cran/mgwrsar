#' bisquare kernel
#' @usage bisq(d, h)
#' @param d  a vector of distance
#' @param h  a distance bandwidth
#' @return a vector of weight
#' @examples
#' \donttest{
#' w=bisq(-30:30, 20)
#' plot(-30:30,w,type='l')
#' abline(v=-20)
#' abline(v=20)
#' }
bisq=function (d, h)
{
  x <- d/h
  x[abs(x) < 1] <- (15/16) * ((1 - x[abs(x) < 1]^2)^2)
  x[abs(x) >= 1] <- 0
  x
}
