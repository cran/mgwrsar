#' kernel_matW
#' A function that returns a sparse weight matrix based computed with a specified
#' kernel (gauss,bisq,tcub,epane,rectangle,triangle) considering coordinates
#' provides in S and a given bandwidth. If NN<nrow(S) only NN firts neighbours are considered.
#' If Type!='GD' then S should have additional columns and several
#' kernels and bandwidths should be be specified by the user.
#' @usage kernel_matW(H,kernels,coords,NN,TP=NULL,Type='GD',adaptive=FALSE,
#' diagnull=TRUE,alpha=1,theta=1,dists=NULL,indexG=NULL,extrapol=FALSE,QP=NULL,K=0)
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param coords A matrix with  variables used in kernel (reference)
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP A vector with index of target points
#' @param Type  Type of Genelarized kernel product ('GD' only spatial,'GDC'
#' spatial + a categorical variable,'GDX' spatial + a continuous variable,
#' 'GDT' spatial + a time index, and other combinations 'GDXXC','GDTX',...)
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel
#' @param diagnull  Zero on diagonal, default FALSE
#' @param  alpha TO BE DOCUMENTED
#' @param theta TO BE DOCUMENTED
#' @param dists TO BE DOCUMENTED
#' @param indexG TO BE DOCUMENTED
#' @param extrapol TO BE DOCUMENTED
#' @param QP  A matrix with  variables used in kernel (neighbors), default NULL (if NULL coord_j=coord_i)
#' @param K TO BE DOCUMENTED
#' @return A sparse Matrix of weights (dgCMatrix).
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix) of 4 nearest neighbors with 0 in diagonal
#'  W=kernel_matW(H=4,kernels='rectangle',coords=coords,NN=4,adaptive=TRUE,diagnull=TRUE)
#' }
kernel_matW<-function(H,kernels,coords,NN,TP=NULL,Type='GD',adaptive=FALSE,diagnull=TRUE,alpha=1,theta=1,dists=NULL,indexG=NULL,extrapol=FALSE,QP=NULL,K=0){
  if(Type=='GDT'){
    if(kernels[1]!='gauss' | unlist(str_split(kernels[2], '_'))[1]!='gauss' ) stop('For Type=GDT only gauss kernel should be used')
  }
 n<-nrow(coords)
  if(is.null(TP)) TP=1:n
  ntp=length(TP)
  if(extrapol) m=n-ntp else m=n
 while(sum(duplicated(coords[,1:2]))>0) {
   coords[,1:2]<-jitter(coords[,1:2],0.001)
      #warning('coords have been jittered because there is some duplicated location.')
    }
  #if(adaptive[1] & Type=='GD' & diagnull==TRUE) NN=min(max(max(H)+1,NN),ntp) ## case locally varying bandwidth
  if(adaptive[1] & Type!='GD' & diagnull==TRUE) NN=min(max(H[1]+2,NN),ntp)
  if(kernels[1]=='shepard') NN=max(H[1]+2,min(NN+1,ntp))
  for(j in 1:length(adaptive)) {
    if(adaptive[j]) NN=max(H[1]+2,min(NN,ntp)) # min(max(NN,H[j]+1),ntp)
    #cat('\n NN increased to ', H[j]+2,' because H[',j,']=',H[j],sep='')
  }
  if((ncol(coords)>2 | Type!='GD') & (length(kernels)<ncol(coords)-2 | length(H)<ncol(coords)-2)) stop("if Type!='GD', H and kernels must be indicated for each variables used in General Kernel Product, i.e length(H)=length(kernels)=ncol(S)-2")
  if((ncol(coords)>2 | Type!='GD') & (length(adaptive)==1)) if(adaptive) stop("if Type!='GD',adaptive must be indicated for each variables used in General Kernel Product, i.e length(adaptive)=ncol(S)-2")

  #if(Type=='GDT' & length(TP)<n) QP=(1:n)[-TP]
  if(is.null(indexG)) stage0=prep_d(coords,NN=NN,TP=TP,extrapol=extrapol,ratio=1,QP=QP)
  indexG=stage0$indexG
  dists=stage0$dists
  stage1=prep_w(H=H,kernels=kernels,Type=Type,adaptive=adaptive,dists=dists,indexG=indexG,alpha=alpha,theta=theta,K=K)
W <- sparseMatrix(i = rep(1:nrow(stage1$indexG),each=ncol(stage1$indexG)), j = t(stage1$indexG),  dims = c(nrow(stage1$indexG),m), x =as.numeric(t(stage1$Wd)))

  if(diagnull) {
    diag(W)=0
    W <- normW(W)
  }
  return(W)
}
