#' kdist_adapt_sorted
#' to be documented
#' @usage kdist(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
kdist_adapt_sorted<-function (d, h) {
  if(is.unsorted(d[1,])) {
    dd=t(apply(d,1,sort))
  } else dd=d
  dd
}
