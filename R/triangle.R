#' triangle
#' to be documented
#' @usage triangle(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
triangle <-function(d,h) {              #triangle function
  x<-abs(d)/h
  x[x<= 1]<-(1-x[x<=1])
  x[x> 1]<- 0
  x
}
