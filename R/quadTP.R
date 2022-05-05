#' quadTP
#' to be documented
#' @usage quadTP(xy,k=100)
#' @param xy to be documented
#' @param k to be documented
#' @noRd
#' @return to be documented
quadTP<-function(xy,k=100){
  xy<-as.matrix(xy)
  colnames(xy)<-c('x','y')
  qt <- quadtree(as.matrix(xy), k)
  xylim <- cbind(x=c(min(xy[,1]), max(xy[,1])), y=c(min(xy[,2]), max(xy[,2])))
  polys <- cell(qt, xylim)
  polys.attr <- data.frame(id=unique(polys$id))
  polys$id<-as.numeric(factor(polys$id))
  insidecell(polys,as.matrix(xy))
}
