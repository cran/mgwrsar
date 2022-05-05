#' equidistantGrid
#' to be documented
#' @usage equidistantGrid(nTP,coord)
#' @param nTP to be documented
#' @param coord to be documented
#' @noRd
#' @return to be documented
equidistantGrid<-function(nTP,coord){
  xmin=min(coord[,1])
  xmax=max(coord[,1])
  ymin=min(coord[,2])
  ymax=max(coord[,2])
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
  index<-unique((knn(coord,coordg,1))$nn.idx)
  #if(length(index)>nTP) index <- index[sample(1:length(index),nTP)]
  sort(index)
}
