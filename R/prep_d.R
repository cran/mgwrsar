#' prep_d
#' to be documented
#' @usage prep_d(coords,NN,TP,extrapol=FALSE,ratio=1,QP=NULL)
#' @param coords A matrix with spatial coordinates
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP index of target points, default 1:n
#' @param extrapol  special mode for prediction unig GWR estimation, default FALSE
#' @param QP index of query points, default NULL
#' @param ratio numeric [0,1] ratio time/space for ordering indexG, default 1.
#' @noRd
#' @return A list with precomputed matrices of distances and neighbours indexes
prep_d<-function(coords,NN,TP,extrapol=FALSE,ratio=1,QP=NULL){
  dists=list()
  if(is.null(QP)) QP=(1:nrow(coords))
  TPc=(1:nrow(coords))[-TP]
  if(ncol(coords)>2){
    #range01 <- function(x){(x-min(x))/(max(x)-min(x))};
    #ki=min(NN,nrow(coords))
    #coords2=apply(coords,2,range01)
    #coords2[,3]=coords2[,3]*ratio
    nn=knn(coords[QP,1:2],k=NN,query=coords[TP,1:2])
    indexG=nn$nn.idx
    dt=matrix(coords[indexG,3],ncol=ncol(indexG),nrow=nrow(indexG))
    dtr=matrix(coords[TP,3],ncol=ncol(indexG),nrow=nrow(indexG),byrow=F)
    dists[['dist_s']]=nn$nn.dists
    dists[['dist_t']]=dtr-dt
  } else {
    if(extrapol) nn=knn(coords[-TP,],k=min(c(NN,nrow(coords)-length(TP))),query=coords[TP,]) else nn=knn(coords,k=NN,query=coords[TP,])
    indexG=nn$nn.idx #    indexG=matrix(TPc[nn$nn.idx],nrow=nrow(nn$nn.idx),ncol=ncol(nn$nn.idx))
    dists[['dist_s']]=nn$nn.dists
  }
  list(indexG=indexG,dists=dists)
}

