#' gauss
#' to be documented
#' @usage gauss(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
gauss<-function (d, h) {
  exp(-0.5*(d/h)^2)
}
