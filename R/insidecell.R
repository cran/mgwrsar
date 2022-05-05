#' insidecell
#' to be documented
#' @usage insidecell(poly,xy)
#' @param poly to be documented
#' @param xy to be documented
#' @noRd
#' @return to be documented
insidecell<-function(poly,xy){
  xy<-data.frame(xy)
  xy$id=0
  for(i in unique(poly$id)){
    pol<-poly[poly$id==i,]
    xy[xy$x>=min(pol$x) & xy$x<=max(pol$x) & xy$y>=min(pol$y) & xy$y<=max(pol$y),'id']<-i
  }
  xy
}
