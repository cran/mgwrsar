#' equidistantGrid
#' to be documented
#' @usage equidistantGrid(nTP,coords)
#' @param nTP to be documented
#' @param coords to be documented
#' @noRd
#' @return to be documented
equidistantGrid<-function(nTP,coords){
  xmin=min(coords[,1])
  xmax=max(coords[,1])
  ymin=min(coords[,2])
  ymax=max(coords[,2])
  rx=xmax-xmin
  ry=ymax-ymin
  surface=(rx*ry)
  pas=sqrt(surface/(nTP-4))
  nx=round(rx/pas)+1
  ny=round(ry/pas)+1
  xx=xmin+((1:nx)-0.5)*pas
  yy=xmin+((1:ny)-0.5)*pas
  coordg=expand.grid(xx,yy)
  colnames(coordg)=c('x','y')
  #coordg<-coordg[sort(sample(1:nrow(coordg),nTP)),]
  index<-unique((knn(coords,coordg,1))$nn.idx)
  #if(length(index)>nTP) index <- index[sample(1:length(index),nTP)]
  sort(index)
}
