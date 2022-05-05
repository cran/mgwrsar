#' bisq
#' to be documented
#' @usage bisq(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
bisq <- function(d,h) {     #bisquare kernel
  x<-d/h
  x[abs(x)< 1]<- (15/16)*((1-x[x< 1]^2)^2)
  x[abs(x)>= 1]<- 0
  x
}
