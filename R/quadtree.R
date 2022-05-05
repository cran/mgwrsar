#' quadtree
#' to be documented
#' @usage quadtree(xy, k=1)
#' @param xy to be documented
#' @param k to be documented
#' @noRd
#' @return to be documented
quadtree<-function(xy, k=1) {
  d = dim(xy)[2]
  quad <- function(xy, i, id=1) {
    if (length(xy) < 2*k*d) {
      rv = list(id=id, value=xy)
      class(rv) <- "quadtree.leaf"
    }
    else {
      q0 <- (1 + runif(1,min=-1/2,max=1/2)/dim(xy)[1])/2 # Random quantile near the median
      x0 <- quantile(xy[,i], q0)
      j <- i %% d + 1 # (Works for octrees, too...)
      rv <- list(index=i, threshold=x0,
                 lower=quad(xy[xy[,i] <= x0, ], j, id*2),
                 upper=quad(xy[xy[,i] > x0, ], j, id*2+1))
      class(rv) <- "quadtree"
    }
    return(rv)
  }
  quad(xy, 1)
}
