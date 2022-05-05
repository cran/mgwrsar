#' prep_d
#' to be documented
#' @usage prep_d(coord,NN,TP)
#' @param coord A matrix with spatial coordinates
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP index of target points, default 1:n
#' @noRd
#' @return A list with precomputed matrices of distances and neighbours indexes
prep_d<-function(coord,NN,TP){
  ### distance and rank matrices
  nn=knn(coord,k=NN,query=coord[TP,])
  indexG=nn$nn.idx
  dists=list(coord=nn$nn.dists)
  list(indexG=indexG,dists=dists$coord)
}
