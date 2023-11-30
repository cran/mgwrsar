#' multiscale_gwr
#' This function adapts the multiscale Geographically Weighted Regression (GWR) methodology
#' proposed by Fotheringam et al. in 2017, employing a backward fitting procedure within
#' the MGWRSAR subroutines. The consecutive bandwidth optimizations are performed by
#' minimizing the corrected Akaike criteria.
#' @usage multiscale_gwr(formula,data,coords,Model = 'GWR',kernels='bisq',
#' control=list(SE=FALSE,adaptive=TRUE,NN=800,isgcv=FALSE),init='GWR',maxiter=100,
#' nstable=6,crit=0.000001,doMC=FALSE,ncore=1,HF=NULL,H0=NULL,model=NULL)
#' @param formula A formula.
#' @param data A dataframe.
#' @param coords default NULL, a dataframe or a matrix with coordinates.
#' @param Model The type of model: Possible values are "GWR" (default),
#' and "MGWRSAR_1_0_kv". See Details for more explanation.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"),
#' gaussian("gauss")).
#' @param control a list of extra control arguments, see MGWRSAR help.
#' @param init starting model (lm or GWR)
#' @param maxiter maximum number of iterations in the back-fitting procedure.
#' @param nstable required number of consecutive unchanged optimal bandwidth (by covariate)
#' before leaving optimisation of bandwidth size, default 3.
#' @param crit value to terminate the back-fitting iterations (ratio of change in RMSE)
#' @param doMC A boolean for Parallel computation, default FALSE.
#' @param ncore number of CPU cores for parallel computation, default 1.
#' @param HF if available, a vector containing the optimal bandwidth parameters for each
#' covariate, default NULL.
#' @param H0 A bandwidth value for the starting GWR model, default NULL.
#' @param model A previous model estimated using multiscale_gwr function, default NULL.
#' @return  Return an object of class mgwrsar with at least the following components:
#' \describe{
#' \item{Betav}{ matrix of coefficients of dim(n,kv) x kv.}
#' \item{Betac}{ vector of coefficients of length kc.}
#' \item{Model}{ The sum of square residuals.}
#' \item{Y}{ The dependent variable.}
#' \item{XC}{ The explanatory variables with constant coefficients.}
#' \item{XV}{ The explanatory variables with varying coefficients.}
#' \item{X}{ The explanatory variables.}
#' \item{W}{ The spatial weight matrix for spatial dependence.}
#' \item{isgcv}{ if gcv has been computed.}
#' \item{edf}{ The estimated degrees of freedom.}
#' \item{formula}{The formula.}
#' \item{data}{ The dataframe used for computation.}
#' \item{Method}{ The type of model.}
#' \item{coords}{ The spatial coordinates of observations.}
#' \item{H}{ A vector of bandwidths.}
#' \item{fixed_vars}{ The names of constant coefficients.}
#' \item{kernels}{ The kernel vector.}
#' \item{SSR}{ The sum of square residuals.}
#' \item{residuals}{ The vector of residuals.}
#' \item{fit}{ the vector of fitted values.}
#' \item{sev}{ local standard error of parameters.}
#' \item{get_ts}{Boolean, if trace of hat matrix Tr(S) should be stored.}
#' \item{NN}{ Maximum number of neighbors for weights computation}
#' }
#' @seealso  tds_mgwr, bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar
#' @examples
#' \donttest{
#' library(mgwrsar)
#' mysimu<-simu_multiscale(n=1000)
#' mydata=mysimu$mydata
#' coords=mysimu$coords
#' model_multiscale<-multiscale_gwr(formula=as.formula('Y~X1+X2+X3'),data=mydata,
#' coords=coords,Model = 'GWR',kernels='bisq',control=list(SE=FALSE,
#' adaptive=TRUE,NN=900,isgcv=FALSE),init='GWR',nstable=6,crit=0.000001)
#' summary_mgwrsar(model_multiscale)
#' }
multiscale_gwr<-function(formula,data,coords,Model = 'GWR',kernels='bisq',control=list(SE=FALSE,adaptive=TRUE,NN=800,isgcv=FALSE),init='GWR',maxiter=100,nstable=6,crit=0.000001,doMC=FALSE,ncore=1,HF=NULL,H0=NULL,model=NULL){
  start<-proc.time()
  n=nrow(data)
  mf <- model.frame(formula, data)
  mt <- attr(x = mf, which = "terms")
  Y <- model.extract(mf, "response")
  X = model.matrix(object = mt, data = mf)
  namesX =colnames(X)
  if(colnames(X)[1]== "(Intercept)") colnames(X)[1]<-namesX[1]<-'Intercept'
  data$Intercept=rep(1,n)
  K=length(namesX)

  if(!is.null(model)){
    HF=model$H
    BETA=model$Betav
    data$eps<-model$residuals
  } else if(init=='lm'){
    model0=lm(formula,data)
    BETA=matrix(rep(coef(model0),each=nrow(data)),byrow=FALSE,ncol=length(coef(model0)))
    data$eps<-residuals(model0)
    H=rep(n,K)
  } else {
    if(is.null(H0)) {
      res=golden_search_bandwidth_AICc(formula = formula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,Model=Model,control=control,lower.bound=5, upper.bound=control$NN-2)
     H0= res$minimum
    }
    modelGWR<-MGWRSAR(formula = formula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=H0, Model = Model,control=control)
    BETA=modelGWR$Betav
    H=rep(modelGWR$H,K)
    data$eps<-modelGWR$residuals
    cat('H0=',H0,'\n')
  }

  drmse<-delta_rmse<-sqrt(mean(Y^2))
  isgcv=FALSE
  iter=0
  G<-prep_d(coords,control$NN,TP=1:nrow(coords))
  control$dists=G$dists
  control$indexG=G$indexG
  stable<-rep(1,K)
  if(!is.null(HF)) H=HF

  while((abs(delta_rmse)>crit | any(stable<nstable)) & iter<=maxiter){
    iter=iter+1
    cat('\n')
    for(k in 1:K){
      var = namesX[k]
      cat(' ', var,' ')
      data$epst<-data$eps+BETA[,k]*data[,var] %>% as.matrix(ncol=1)
      myformula=as.formula(paste0('epst~',var,'-1'))
      if(is.null(HF)){
        if(stable[k]<nstable){
          res=golden_search_bandwidth_AICc(formula = myformula, data = data,coords=coords, kernels=kernels,fixed_vars=NULL,Model=Model,control=control,lower.bound=5, upper.bound=control$NN-2)
          modellm=lm(formula = myformula, data = data)
          if(res$objective>AIC(modellm)){
            res$minimum=n
            res$objective=AIC(modellm)
          }
          H[k]<-h<-res$minimum } else h=H[k]} else h=HF[k]

      if(h==n){
        modelGWR<-lm(myformula,data)
        modelGWR$H=n
        modelGWR$Betav=rep(coef(modelGWR),n)
      } else {
        if(stable[k]<nstable & is.null(HF)) modelGWR<-res$model else modelGWR<-MGWRSAR(formula =myformula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=h, Model = Model,control=control)
      }
      BETA[,k]<-modelGWR$Betav
      data$eps=modelGWR$residuals
      if(iter>1) if(oldH[k]==H[k]) stable[k]=stable[k]+1 else stable[k]=1
    }
    oldH<-H
    delta_rmse=(drmse - sqrt(mean(data$eps^2)))/drmse
    drmse= sqrt(mean(data$eps^2))
    cat('\n iter ',iter,' H ',H,' RMSE=',sqrt(mean(data$eps^2)),' delta_rmse ',delta_rmse,' stable',stable,'\n')
  }
  cat('Duree = ',(proc.time()- start)[3],'\n')
  modelGWR$Model='multiscale_gwr'
  modelGWR$Betav=BETA
  modelGWR$residuals=data$eps
  modelGWR$fitted=Y-data$eps
  modelGWR$RMSE=sqrt(mean(modelGWR$residuals^2))
  modelGWR$H=H
  modelGWR$X=X
  modelGWR
}
