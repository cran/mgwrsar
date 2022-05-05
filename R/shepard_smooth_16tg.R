#' An extrapolation method using shepard smooth with a given number of neighbors.
#' @usage Beta_extropolation(coord=coord,k=16,O=NULL,type='SxS')
#' @param coord a matrix of coordinates
#' @param k The number of neighbors for extrapolation
#' @param O An index vector for outsample
#' @param type Insample/outsample size configuration : 'S+OxS+O' or 'SxS'(default)
#' @noRd
#' @return A weight dgCMatrix of size $type$
Beta_extropolation <-
function(coord=coord,k=16,O=NULL,type='SxS'){
if(is.null(O) & type=='S+OxS+O') stop('You have to priovide O when type= "S+OxS+O"')
n=nrow(coord)
S=1:n
m=n
if(!is.null(O)) {S=S[-O]
m=length(O)
type='S+OxS+O'
}
if(length(S)<k) stop(paste('choose at least a insample size >=', k))

if (type=='SxS') {
nb1<-knn(as.matrix(coord), query =  as.matrix(coord), k = k+1)
nb1$nn.dists<-nb1$nn.dists[,-1]
nb1$nn.idx<-nb1$nn.idx[,-1]
} else {
nb1<-knn(as.matrix(coord)[S,], query =  as.matrix(coord)[O,], k = k)
}
DMmax<-apply(nb1$nn.dists,1,max)
d<-as.numeric(t(nb1$nn.dists))
dmax<-as.numeric(sapply(DMmax,function(x) rep(x,k)))
dd=(d-dmax)^2/(d*dmax)
if (type=='SxS'){
W<-sparseMatrix(i=rep(1:n,each=k),j=t(nb1$nn.idx),x=dd)
W<-normW(W)
} else if(type=='S+OxS+O')
{
index<-matrix(S[nb1$nn.idx],nrow = nrow(nb1$nn.idx), ncol = ncol(nb1$nn.idx),byrow = FALSE)
W<-sparseMatrix(i=O[rep(1:m,each=k)],j=t(index),x=dd ,dims=c(n,n),symmetric=FALSE)
W<-normW(W)
} # else if (type=='S+OxS')
# {
# W<-sparseMatrix(i=O[rep(1:m,each=k)],j=t(nb1$nn.idx),x=dd ,dims=c(n,n-m),symmetric=FALSE)
# rs <- W %*% as.matrix(rep(1, n-m))
# W[rs>0]<-W[rs>0]/rs[rs>0]
# } else {
# W<-sparseMatrix(i=rep(seq_along(kk),kk),j=t(nb1$nn.idx),x=dd ,dims=c(n-m,m),symmetric=FALSE)
# rs <- W %*% as.matrix(rep(1, m))
# W[rs>0]<-W[rs>0]/rs[rs>0]
# }
##W<-as.matrix(W)
if(m!=n) Matrix::diag(W)[S]=1
W
}
