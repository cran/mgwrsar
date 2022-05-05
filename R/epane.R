#' epane
#' to be documented
#' @usage epane(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
epane <-function(d,h) {              #epanechnikov kernel
  x<-abs(d)/h
  x[x< 1]<- 3/4*(1-(x[x< 1])^2)
  x[x>= 1]<- 0
  x
}
