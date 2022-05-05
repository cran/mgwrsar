#' prep_w
#' to be documented
#' @usage prep_w(H,kernels,coord_i,coord_j,NN,ncolX,Type='GD',adaptive=F,dists=NULL,
#' indexG=NULL,rowNorm=TRUE,extrapTP=0,correctionNB=TRUE)
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param coord_i  A matrix with spatial coordinates and eventually other variables  (ref)
#' @param coord_j  A matrix with spatial coordinates and eventually other variables (neighbors), default NULL.
#' @param NN Number of spatial Neighbours for kernels computations
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel
#' @param ncolX  The number of model covariates
#' @param dists  Precomputed Matrix of spatial distances, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL
#' @param rowNorm  A boolean, row normalization of weights, default TRUE
#' @noRd
#' @return to be documented
prep_w<-function(H,kernels,coord_i,coord_j,NN,ncolX,Type='GD',adaptive=F,dists=NULL,indexG=NULL,rowNorm=TRUE){
  if(ncol(coord_i)!=ncol(coord_j)) stop("coord_i and coord_j must have the same number of columns")
  n=nrow(coord_j)
  ntp=nrow(coord_i)
  if(ncol(coord_i)>2) {
    Z=as.matrix(coord_i[,-(1:2)])
    Z_in=as.matrix(coord_j[,-(1:2)])
    coord_i<-as.matrix(coord_i[,(1:2)])
    coord_j<-as.matrix(coord_j[,(1:2)])
    } else Z=NULL

  ### adpative : if adaptive round(H) and rename kernel
  for(j in 1:length(adaptive)){
    case = adaptive[j]
    if(case) {
      if(Type=='GD' & any(H!=H[1])) H=round(H) else H[j]=round(H[j])
      kernels[j]= paste0(kernels[j],'_adapt_sorted')
    }
  }

  myH=H
  if(is.null(colnames(Z)) & !is.null(Z)) colnames(Z)=paste('Z',1:ncol(Z))

  if(length(H)<ntp) {
    H=list(coord=rep(H[1],ntp)) # ifelse(extrapTP==0,n,ntp)
    if(nchar(Type)>2) for(t in 1:(nchar(Type)-2)) {
      if(substr(Type,t+2,t+2)=='C'){
        nbcs=max(Z[,t])
        for(nbc in 1:nbcs){
          H[[paste0(colnames(Z)[t],'_',nbc)]]=rep(myH[t+nbc],ntp)
        }
      } else  H[[colnames(Z)[t]]]=rep(myH[t+1],ntp)
    }
  } else H=list(coord=H)
  ### distance and rank matrices
  if(is.null(dists)){

    # if(extrapTP==0) { # here we want length(TP) x nrow(X) weight matrix for estimation on Betav(TP)
    #   nn=knn(coord,k=NN,query=coord[TP,])
    #   indexG=nn$nn.idx
    # } else if(extrapTP==1) {# here we want nrow(X) x length(TP) weight matrix for extrapolating Betav(-TP) using Betav(TP), this function return a n x n matrix and the selection of -TP lines is done after, outside this function.
    #   ## option 1 : projection classique
    #   nn=knn(coord[TP,],k=NN,query=coord)
    #   indexG=matrix(TP[nn$nn.idx],ncol=ncol(nn$nn.idx),nrow=nrow(nn$nn.idx))
    # } else if(extrapTP==2) {# here we want nrow(rbind(S,O)) x length(S) weight matrix for extrapolating Betav(O) using Betav(S), this function return a n x n matrix and the selection of -TP lines is done after, outside this function.
    #   ## option 1 : projection classique
    #   nn=knn(coord[TP,],k=NN,query=coord)
    #   indexG=matrix(TP[nn$nn.idx],ncol=ncol(nn$nn.idx),nrow=nrow(nn$nn.idx))
    # } else if(extrapTP==3){
    #   nn=knn(coord[TP,],k=min(NN,length(TP)-1),query=coord[-TP,])
    #   indexG=nn$nn.idx
    # }

    ####
    ## if no TP estimation : coord_i = coord_j=coord , matrix W n x n
    ## if TP estimation : coord_i=coord[TP] \in coord_j=coord , matrix W ntp x n
    ## if m extrapolation from Beta : coord_i \not in coord_j=new_data, matrix W n x m
    ## if m extrapolation from Beta(TP) : coord_i=coord[TP] \not in coord_j=new_data, matrix W ntp x m

    nn=knn(coord_j,k=min(NN,n),query=coord_i)
    indexG=nn$nn.idx
    dists=list(coord=nn$nn.dists)
  } else dists=list(coord=dists)
  ## GPK with D
  mykernels=kernels
  kernels=list(coord=mykernels[1])
  Wd=do.call(kernels$coord,args=list(dists$coord,H$coord))
  if(rowNorm) Wd=Wd/rowSums(Wd)
   nv=apply(Wd,1,function(x) sum(x>0))
   ### case non adaptive with locally not enough neighbors :
   if(any(nv<2*ncolX & !adaptive[1] & kernels$coord!='sheppard' )){
     index=which(nv<2*ncolX)
     Wd[index,]<-do.call(paste0(kernels$coord,'_adapt_sorted'),args=list(matrix(dists$coord[index,],ncol=ncol(dists$coord)),rep(2*ncolX,length(index))))
   }
  if(nchar(Type)>2){
    cases='D'
    for(t in 1:(nchar(Type)-2)) {
      kernels[[colnames(Z)[t]]]=mykernels[t+1]
      typek=substr(Type,t+2,t+2)
      if(typek!='C'){
      if(sum(duplicated(Z[,t]))>0) Z[,t]=jitter(Z[,t],0.001)
      dists[[colnames(Z)[t]]]=t(apply(cbind(1:nrow(indexG),indexG),1, function(x) abs(Z[x[1],t]-Z_in[x[-1],t]) ))
      wd=do.call(kernels[[colnames(Z)[t]]],args=list(dists[[colnames(Z)[t]]],H[[colnames(Z)[t]]]))
      } else {
        wd=matrix(NA,nrow=ntp,ncol=NN)
        for(i in 1:ntp){
          wd[i,] <- ifelse(Z[indexG[i,],t] != Z[i,t], myH[[as.numeric(t + Z[i,t])]], 1)
        }

      }
      # if(typek=='T'){
      #     post<-t(apply(indexG,1, function(x) {
      #       post<-as.numeric(Z[x[1],t]>=Z[x,t])
      #       post[post==0]<-0.3
      #       post
      #       }))
      #     wd<-wd*post
      # }
      if(rowNorm) wd=wd/rowSums(wd)
      cases=c(cases,typek)
      Wd=Wd*wd
    }
  }
  if(rowNorm) Wd=Wd/rowSums(Wd)
  list(indexG=indexG,Wd=Wd,dists=dists$coord)
}
