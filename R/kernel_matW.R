#' kernel_matW
#' A function that returns a sparse weight matrix based computed with a specified
#' kernel (gauss,bisq,tcub,epane,rectangle,triangle) considering coordinates
#' provides in S and a given bandwidth. If NN<nrow(S) only NN firts neighbours are considered.
#' If Type!='GD' then S should have additional columns and several
#' kernels and bandwidths should be be specified by the user.
#' @usage kernel_matW(H,kernels,coord_i,coord_j=NULL,NN,ncolX=1,
#' Type='GD',adaptive=F,diagnull=TRUE,rowNorm=TRUE,noisland=FALSE)
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param coord_i  A matrix with  variables used in kernel (reference)
#' @param coord_j  A matrix with  variables used in kernel (neighbors), default NULL (if NULL coord_j=coord_i)
#' @param NN Number of spatial Neighbours for kernels computations
#' @param ncolX control parameter
#' @param Type  Type of Genelarized kernel product ('GD' only spatial,'GDC'
#' spatial + a categorical variable,'GDX' spatial + a continuous variable,
#' 'GDT' spatial + a time index, and other combinations 'GDXXC','GDTX',...)
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel
#' @param diagnull  Zero on diagonal, default FALSE
#' @param rowNorm  A boolean, row normalization of weights, default TRUE
#' @param noisland A boolean to avoid isle with no neighbours for non adaptive kernel, default FALSE
#' @return A sparse Matrix of weights (dgCMatrix).
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coord=as.matrix(mydata[,c("x_lat","y_lon")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix) of 4 nearest neighbors with 0 in diagonal
#'  W=kernel_matW(H=4,kernels='rectangle',coord_i=coord,NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=TRUE)
#' }
kernel_matW<-function(H,kernels,coord_i,coord_j=NULL,NN,ncolX=1,Type='GD',adaptive=F,diagnull=TRUE,rowNorm=TRUE,noisland=FALSE){
  m<-n<-nrow(coord_i)
  if(is.null(coord_j)){
    while(sum(duplicated(coord_i))>0) {
      coord_i<-jitter(coord_i,0.001)
      warning('coords have been jittered because there is some duplicated location.')
    }
    coord_j=coord_i
  } else m<-nrow(coord_j)
  if(ncol(coord_i)!=ncol(coord_j)) stop("coord_i and coord_j must have the same number of columns")
  if(adaptive[1] & Type=='GD' & diagnull==TRUE) NN=min(max(max(H)+1,NN),m) ## case locally varying bandwidth
  if(adaptive[1] & Type!='GD' & diagnull==TRUE) NN=min(max(H[1]+1,NN),m)
  if(kernels[1]=='sheppard') NN=min(NN+1,m)
  for(j in 1:length(adaptive)) {
    if(adaptive[j]) NN=min(max(NN,H[j]+1),m)
    #cat('\n NN increased to ', H[j]+2,' because H[',j,']=',H[j],sep='')
  }
  if((ncol(coord_i)>2 | Type!='GD') & (length(kernels)<ncol(coord_i)-2 | length(H)<ncol(coord_i)-2)) stop("if Type!='GD', H and kernels must be indicated for each variables used in General Kernel Product, i.e length(H)=length(kernels)=ncol(S)-2")
  if((ncol(coord_i)>2 | Type!='GD') & (length(adaptive)==1)) if(adaptive) stop("if Type!='GD',adaptive must be indicated for each variables used in General Kernel Product, i.e length(adaptive)=ncol(S)-2")

  stage1=prep_w(H=H,kernels=kernels,coord_i=coord_i,coord_j=coord_j,NN=NN,ncolX=ncolX,adaptive=adaptive,dists=NULL,indexG=NULL,rowNorm=TRUE,noisland=noisland)
W <- sparseMatrix(i = rep(1:nrow(stage1$indexG),each=ncol(stage1$indexG)), j = t(stage1$indexG),  dims = c(nrow(stage1$indexG),m), x =as.numeric(t(stage1$Wd)))

  if(diagnull) {
    diag(W)=0
    if(rowNorm) W <- normW(W)
  }
  return(W)
}
