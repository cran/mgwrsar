#' Search of a suitable set of target points.
#' tp is a function that identifies a set of target points based on spatial smoothed OLS residuals.
#' @usage tp(K,Residuals,coord,kWtp,Wtp=NULL,prev_TP=NULL)
#' @param K the minimum number of first neighbors with lower (resp.higer) absolute value of the smoothed residuals.
#' @param Residuals a vector of residuals
#' @param coord a dataframe or a matrix with coordinates, not required if data is a spatial dataframe
#' @param kWtp the number of first neighbors for computing  the smoothed residuals, default 16.
#' @param prev_TP index of already used TP (version length(K)>1), default NULL.
#' @return a list with two vectors, the index of target points and a vector of corresponding smoothed residuals.
#' @details tp is a function that identifies a set of target points based on spatial smoothed residuals.
#' The function first computes the smooth of model residuals using a Sheppard's kernel with kWtp neighbors (default 16).
#' Then it identifies local maxima (resp. minima) that fits the requirement of having at least K neighbors with lower (resp.higer) absolute value of the smoothed residuals. As K increases the number of target points decreases.
#' @noRd
tp<-function(K,Residuals,coord,kWtp,Wtp=NULL,prev_TP=NULL){
  n=length(Residuals)
  indexG<-(knn(coord,k=max(kWtp+1,K+1)))$nn.idx
  Mr<-matrix(Residuals[indexG[,2:(kWtp+1)]],ncol=kWtp) ## i js residuals
  #Wres=apply(Mr,1,sum)/kWtp
  Wres=rowMeans(Mr)  ## row mean of i js residuals
  Mwr<-matrix(abs(Wres)[indexG[,1:K]],ncol=K)   ##
  if(K!=n){
    tgmaxmin=apply(Mwr,1,which.max) ##tgmaxmin=Rfast::rowMaxs(Mwr,F)
    tgmaxmin=which(tgmaxmin==1)
  if(!is.null(prev_TP)) tgmaxmin=setdiff(tgmaxmin,prev_TP)} else tgmaxmin=which.max(abs(Wres))
  list(tgmaxmin=tgmaxmin,vmaxmin= Wres[tgmaxmin])
}
