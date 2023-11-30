#' prep_d
#' to be documented
#' @usage prep_d(coords,NN,TP)
#' @param coords A matrix with spatial coordinates
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP index of target points, default 1:n
#' @noRd
#' @return A list with precomputed matrices of distances and neighbours indexes
prep_d<-function(coords,NN,TP){
  ### distance and rank matrices
  nn=knn(coords,k=NN,query=coords[TP,])
  indexG=nn$nn.idx
  dists=list(coords=nn$nn.dists)
  list(indexG=indexG,dists=dists$coords)
}
