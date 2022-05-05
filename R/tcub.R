#' tcub
#' to be documented
#' @usage tcub(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @return a vector of weights.
#' @noRd
tcub <- function(d,h) {     #bisquare kernel
  x<-abs(d/h)
  x[x< 1]<- (70/81)*((1-x[x< 1]^3)^3)
  x[x>= 1]<- 0
  x
}
