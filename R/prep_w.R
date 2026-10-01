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

  # temporal_distance_modulo <- function(x, cycling = 365) {
  #   x_mod <- x %% cycling
  #   x_mod[x_mod == 0] <- cycling
  #   pmin(x_mod, cycling - x_mod)
  # }
  # An adaptive bandwidth is a number of neighbours. The compact kernels read
  # the (H + 2)-th sorted distance as bandwidth (rectangle: the (H + 1)-th), so
  # H is capped where that column exists; it used to be capped at the number
  # of columns, one or two short, hence an index error for H close to NN.
  adapt_cap <- function(kernel, ncol_d)
    ncol_d - switch(sub("_.*$", "", kernel), bisq = , epane = , tcub = , triangle = 2L, rectangle = 1L, 0L)
  if(adaptive[1]) {
    H[1]=round(H[1])
    H[1] <- min(H[1], adapt_cap(kernels[1], ncol(indexG)))
    kernels[1]= paste0(kernels[1],'_adapt_sorted')
  }

  # Normalised kernel matrices are memoised per (distance matrix, kernel,
  # bandwidth, format) in the session cache of R/weights_cache.R: identical
  # values, computed once per distinct bandwidth instead of once per call.
  temporal_w <- function(kernels_t, format_t, h, n_norm) {
    .mgwrsar_weights(dists[['dist_t']], list("t", kernels_t, format_t, h, n_norm), function() {
      is_past <- !is.na(format_t) && format_t=='past'
      # the 'past' mask applies between kernel and normalisation
      wt=kernel_eval(kernels_t,dists[['dist_t']],h, if (is_past) 0L else n_norm)
      if(is_past) { ### only past observations for i> H[2]
        past=(dists[['dist_t']]>=0)*1
        id<-which(rowSums(past)>4)
        if(length(id)>0) wt[id,]<-wt[id,]*past[id,] else wt<-wt*past
        for (r in seq_len(n_norm)) wt=normW(wt)
      }
      wt
    })
  }
  spatial_w <- function(n_norm) {
    .mgwrsar_weights(dists[['dist_s']], list("s", kernels[1], H[1], n_norm), function() {
      kernel_eval(kernels[1],dists[['dist_s']],H[1], n_norm)
    })
  }

  if(Type=='T') {
    kernels_t<-unlist(str_split(kernels, '_'))[1]
    # suffix the base kernel name: `kernels` already carries '_adapt_sorted' from
    # the block above (and possibly a '_past' / cycling format suffix)
    if(adaptive[1]) {
      H[1] <- min(round(H[1]), adapt_cap(kernels_t, ncol(dists[['dist_t']])))
      kernels_t= paste0(kernels_t,'_adapt_sorted')
    }
    format_t<-unlist(str_split(kernels, '_'))[2]
    #cycling<-as.numeric(unlist(str_split(kernels, '_'))[3])
    Wd=temporal_w(kernels_t, format_t, H[1], 2L)
  } else if(Type=='GDT') {
    Wd=spatial_w(1L)
    kernels_t<-unlist(str_split(kernels[2], '_'))[1]
    if(adaptive[2]) {
      H[2] <- min(round(H[2]), adapt_cap(kernels_t, ncol(dists[['dist_t']])))
      kernels_t= paste0(kernels_t,'_adapt_sorted')
    }
    format_t<-unlist(str_split(kernels[2], '_'))[2] ## in control ?
    wt=temporal_w(kernels_t, format_t, H[2], 1L)
    if (alpha == 1) {
      # Standard case: Pure product (Interactions), fused with the final normW
      Wd <- .mgwrsar_prod_normW(Wd, wt)
    } else {
      if (alpha == 0) {
        # Rare case: Pure sum (Additive)
        Wd <- Wd + wt
      } else {
        # Mixed case: alpha * (Wd * wt) + (1 - alpha) * (Wd + wt)
        Term_Inter <- Wd * wt
        Wd <- Wd + wt
        Wd <- Wd * (1 - alpha)
        # Final assembly
        Wd <- Wd + (Term_Inter * alpha)
        rm(Term_Inter)
      }
      Wd=normW(Wd)
    }
    rm(wt)
  } else {
    Wd=spatial_w(2L)
  }
  list(indexG=indexG,Wd=Wd,dists=dists)
}
