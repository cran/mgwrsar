#' bandwidth_backfit
#' to be documented (experimental)
#' @usage bandwidth_backfit(formula, data, coord, fixed_vars = NULL, kernels, H,Model = "GWR", isgcv=FALSE,
#' control = list(),fd=c(0.6,0.9,0.99,1,1.01,1.1,1.4))
#' @param formula  a formula.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coord default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param H  A vector of bandwidths
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"), gaussian
#' ("gauss")).
#' @param isgcv  leave one out cross validation, default FALSE
#' @param control a list
#' @noRd
bandwidth_backfit<-function(formula, data, coord, fixed_vars = NULL, kernels, H,Model = "GWR", isgcv=FALSE,control = list(),fd=c(0.6,0.9,0.99,1,1.01,1.1,1.4)){
  n=nrow(data)
  if(control$adaptive){
    myseq1<-1:(which(fd==1)-1)
    myseq2<-(which(fd==1)+1):length(fd)
    Hs<-c(sapply(H*fd[myseq1],floor),H,sapply(H*fd[myseq2],ceiling))
    Hs<-Hs[Hs<n-1]
  } else Hs=H*fd

  control$isgcv=TRUE
  dists=control$dists
  indexG=control$indexG
  verbose=control$verbose
  Hs=unique(Hs)
  res<-foreach(H =Hs,.combine =c,.inorder=FALSE)  %dopar%  {
    model<-MGWRSAR(formula = formula, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=H, Model = Model,control=control)
    model$RMSE
  }
  H=Hs[which.min(res)]
  CV<-res[which.min(res)]
  control$isgcv=isgcv
  model<-MGWRSAR(formula = formula, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=H, Model = Model,control=control)
  if(verbose) cat(Hs,' H=',H,' CV= ',CV,' RMSE= ',model$RMSE,' var=',attr(terms(formula),"term.labels"),'\n')
  list(Betav=model$Betav,fit=model$fit,residuals=model$residuals,CV=CV,RMSE=model$RMSE,H=H)
}
