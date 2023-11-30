#' Search of a suitable set of target points.
#' find_TP is a wrapper function that identifies  a set of target points based on spatial smoothed OLS residuals.
#' @usage find_TP(formula, data,coords,K,kWtp=16,Wtp=NULL,type='residuals',
#' model_residuals=NULL,verbose=0,prev_TP=NULL,nTP=NULL)
#' @param formula a formula
#' @param data a dataframe or a spatial dataframe (SP package)
#' @param coords a dataframe or a matrix with coordinates, not required if data is a spatial dataframe
#' @param K the minimum number of first neighbors with lower (resp.higer) absolute value of the smoothed residuals.
#' @param Wtp a precomputed matrix of weights, default NULL.
#' @param kWtp the number of first neighbors for computing  the smoothed residuals, default 16.
#' @param type method for choosing TP, could be 'residuals', 'equidistantGrid','random',  default 'residuals'
#' @param model_residuals (optional) a vector of residuals.
#' @param verbose verbose mode, default FALSE.
#' @param prev_TP index of already used TP (version length(K)>1), default NULL.
#' @param nTP numbeer of target points for random choice of target points, default NULL.
#' @return  find_TP returns an index vector of Target Points set.
#' @details find_TP is a wrapper function that identifies a set of target points, based on spatial smoothed residuals by default.
#' If no vector of residuals are provided, OLS residuals are computed.
#' The function first computes the smooth of model residuals using a Sheppard's kernel with kWtp neighbors (default 16).
#' Then it identifies local maxima (resp. minima) that fits the requirement of having at least K neighbors with lower (resp.higer) absolute value of the smoothed residuals. As K increases the number of target points decreases.
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  TP=find_TP(formula = 'Y_gwr~X1+X2+X3', data =mydata,coords=coords,K=6,type='residuals')
#'  # only 60 targets points are used
#'  length(TP)
#'
#'  model_GWR_tp<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coords=coords,
#'  fixed_vars=NULL,kernels=c('gauss'),  H=0.03, Model = 'GWR',
#'  control=list(SE=TRUE,TP=TP,kWtp=12))
#'  summary(model_GWR_tp$Betav)
#'  }
find_TP <-function(formula, data,coords,K,kWtp=16,Wtp=NULL,type='residuals',model_residuals=NULL,verbose=0,prev_TP=NULL,nTP=NULL){
  n<-nrow(data)
  if(is.null(nTP)) nTP=round(n/K)
  if(type=='residuals'){
    if(K<=1) stop ("A non-zero positive integer is required for K")
    if(verbose>0) cat('\n-------------------------------------------------\n Search of Target Points \n-------------------------------------------------\n')

    if(K==1) TP=1:n else {
      if(is.null(model_residuals)) model_residuals= residuals(lm(formula,data))
      mycandidats=tp(K,model_residuals,coords,kWtp,Wtp,prev_TP)
      mycandidats=unlist(mycandidats$tgmaxmin)
      TP=as.numeric(mycandidats)
    }
  } else if(type=='equidistantGrid') {
    TP=equidistantGrid(nTP,coords)
  } else if(type=='random') {
    TP=sample(1:n,nTP)
  }
  TP[!duplicated(coords[TP,])]
}

