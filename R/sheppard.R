#' sheppard
#' to be documented
#' @usage sheppard(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of number of neighbours as bandwidths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
sheppard<- function (d, h) {
    if(any(h!=h[1])) {dh <- d[cbind(1:nrow(d),h+1)]} else dh=d[,h[1]+1]
    w=(d-dh)^2/(d*dh)
    ii=which(is.infinite(w[,1]))
    w[ii,1]<-1
    w[ii,2:(ncol(w))]<-0
    jj=which(is.na(w[,1]))
    w[jj,]<-1
    w
  }
