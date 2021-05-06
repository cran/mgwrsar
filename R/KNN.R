#' A function that returns a row normalized weight matrix based on k first neighbors,
#' to be documented
#' @usage KNN(coord,h,diagnull=TRUE,kernel='knn',query=NULL)
#' @param coord  matrix of coordinates
#' @param h  a bandwidth
#' @param diagnull  0 on diagonal, default TRUE
#' @param kernel  kernel type ('bisq','bisq_knn','gauss','gauss_adapt','knn')
#' @param query  an index of neighbors to consider, if NULL all observation are used.
#' @return a row nomralized weight dgCmatrix
#' @examples
#' \donttest{
#' data(mydata)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' W=KNN(coord,8)
#' which(W[1,]>0)
#' W[1,W[1,]>0]
#' }
KNN<-function(coord,h,diagnull=TRUE,kernel='knn',query=NULL){
     if(is.null(query)) query=coord
    if(!is.null(dim(coord))) n=nrow(coord) else n=length(coord)
    h=h+1
    nb1<-knn(as.matrix(coord),query, k = h)
    d <- as.numeric(t(nb1$nn.dists))
    DMmax=nb1$nn.dists[,h]
    DMmax=DMmax*1.01
    dmax <- as.numeric(sapply(DMmax, function(x) rep(x, h))) ## DOMC ?
    dmax[dmax==0]<-1
    if(kernel=='bisq_knn') dd = bisq(d,dmax) else if(kernel=='gauss_adapt') dd = gauss_adapt(d, dmax) else if(kernel=='sheppard') dd=(d - dmax)^2/(d * dmax) else dd=1/h

    m=nrow(query)
    W <- sparseMatrix(i = rep(seq_along(rep(h, m)), rep(h, m)), j = t(nb1$nn.idx),  dims = c(m,n), x = dd)
    if(diagnull) diag(W)=0
    normW(W)
  }
