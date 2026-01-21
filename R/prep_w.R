#' prep_w
#' to be documented
#' @usage prep_w(H,kernels,Type='GD',adaptive=FALSE,dists=NULL,indexG=NULL,alpha=1)
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param Type Type of
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel
#' @param dists  Precomputed Matrix of spatial distances, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL
#' @noRd
#' @return to be documented
prep_w<-function(H,kernels,Type='GD',adaptive=FALSE,dists=NULL,indexG=NULL,alpha=1){

  temporal_distance_modulo <- function(x, cycling = 365) {
    x_mod <- x %% cycling
    x_mod[x_mod == 0] <- cycling
    pmin(x_mod, cycling - x_mod)
  }
  if(adaptive[1]) {
    H[1]=round(H[1])
    kernels[1]= paste0(kernels[1],'_adapt_sorted')
    if(H[1]>ncol(indexG)) H[1]<-ncol(indexG)
  }

  if(Type=='T') {
    kernels_t<-unlist(str_split(kernels, '_'))[1]
    if(adaptive)  kernels_t= paste0(kernels,'_adapt_sorted')
    format_t<-unlist(str_split(kernels, '_'))[2]
    cycling<-as.numeric(unlist(str_split(kernels, '_'))[3])
    wt=do.call(kernels_t,args=list(dists[['dist_t']],H[1]))

    if(!is.na(format_t)){
      if(format_t=='past') { ### only past observations for i> H[2]
        past=(dists[['dist_t']]>=0)*1
        id<-which(rowSums(past)>4)
        if(length(id)>0) wt[id,]<-wt[id,]*past[id,] else wt<-wt*past
      }
    }
    Wd=wt
    } else {
    Wd=do.call(kernels[1],args=list(dists[['dist_s']],H[1]))
  }
  Wd=normW(Wd)
  if(Type=='GDT') {
    kernels_t<-unlist(str_split(kernels[2], '_'))[1]
    if(adaptive[2])  kernels_t= paste0(kernels[2],'_adapt_sorted')
     format_t<-unlist(str_split(kernels[2], '_'))[2] ## in control ?
      wt=do.call(kernels_t,args=list(dists[['dist_t']],H[2]))

    if(!is.na(format_t)){
      if(format_t=='past') {
        past=(dists[['dist_t']]>=0)*1
        id<-which(rowSums(past)>4)
        if(length(id)>0) wt[id,]<-wt[id,]*past[id,] else wt<-wt*past
      }
    }
    wt=normW(wt)
    if (alpha == 1) {
      # Standard case: Pure product (Interactions)
      Wd <- Wd * wt
    } else if (alpha == 0) {
      # Rare case: Pure sum (Additive)
      Wd <- Wd + wt
    } else {
      # Mixed case: alpha * (Wd * wt) + (1 - alpha) * (Wd + wt)
      Term_Inter <- Wd * wt
      Wd <- Wd + wt
      Wd <- Wd * (1 - alpha)
      # Final assembly
      Wd <- Wd + (Term_Inter * alpha)
      # Cleanup
      rm(Term_Inter)
    }
    # Remove wt as soon as possible
    rm(wt)
  }
  Wd=normW(Wd)
  list(indexG=indexG,Wd=Wd,dists=dists)
}
