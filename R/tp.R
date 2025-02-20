#' Search of a suitable set of target points.
#' tp is a function that identifies a set of target points based on spatial smoothed OLS residuals.
#' @usage tp(kt,Residuals,coords,ks,Wtp=NULL,prev_TP=NULL)
#' @param kt the minimum number of first neighbors with lower (resp.higer) absolute value of the smoothed residuals.
#' @param Residuals a vector of residuals
#' @param coords a dataframe or a matrix with coordinates, not required if data is a spatial dataframe
#' @param ks the number of first neighbors for computing  the smoothed residuals, default 16.
#' @param prev_TP index of already used TP (version length(kt)>1), default NULL.
#' @return a list with two vectors, the index of target points and a vector of corresponding smoothed residuals.
#' @details tp is a function that identifies a set of target points based on spatial smoothed residuals.
#' The function first computes the smooth of model residuals using a Shepard's kernel with ks neighbors (default 16).
#' Then it identifies local maxima (resp. minima) that fits the requirement of having at least kt neighbors with lower (resp.higer) absolute value of the smoothed residuals. As kt increases the number of target points decreases.
#' @noRd
tp<-function(kt,Residuals,coords,ks,Wtp=NULL,prev_TP=NULL){
  n=length(Residuals)
  indexG<-(knn(coords,k=max(ks+1,kt+1)))$nn.idx
  Mr<-matrix(Residuals[indexG[,2:(ks+1)]],ncol=ks) ## i js residuals
  #Wres=apply(Mr,1,sum)/ks
  Wres=rowMeans(Mr)  ## row mean of i js residuals
  Mwr<-matrix(abs(Wres)[indexG[,1:kt]],ncol=kt)   ##
  if(kt!=n){
    tgmaxmin=apply(Mwr,1,which.max) ##tgmaxmin=Rfast::rowMaxs(Mwr,FALSE)
    tgmaxmin=which(tgmaxmin==1)
  if(!is.null(prev_TP)) tgmaxmin=setdiff(tgmaxmin,prev_TP)} else tgmaxmin=which.max(abs(Wres))
  list(tgmaxmin=tgmaxmin,vmaxmin= Wres[tgmaxmin])
}
