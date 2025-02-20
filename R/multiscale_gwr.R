#' multiscale_gwr
#' This function adapts the multiscale Geographically Weighted Regression (GWR)
#' methodology proposed by Fotheringam et al. in 2017, employing a backward
#' fitting procedure within the MGWRSAR subroutines. The consecutive bandwidth
#' optimizations are performed by minimizing the corrected Akaike criteria.
#' @usage multiscale_gwr(formula,data,coords,kernels='bisq',init='GWR',
#' maxiter=20,nstable=6,tolerance=0.000001,doMC=FALSE,ncore=1,HF=NULL,
#' H0=NULL,H2=NULL,Model=NULL,model=NULL,get_AICg=FALSE,verbose=FALSE,
#' control=list(SE=FALSE,adaptive=TRUE,NN=800,isgcv=FALSE,family=gaussian()))
#' @param formula A formula.
#' @param data A dataframe.
#' @param coords default NULL, a dataframe or a matrix with coordinates.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane")
#' @param init starting model (lm or GWR)
#' @param maxiter maximum number of iterations in the back-fitting procedure.
#' @param nstable required number of consecutive unchanged optimal bandwidth (by covariate) before leaving optimisation of bandwidth size, default 3.
#' @param tolerance value to terminate the back-fitting iterations (ratio of change in RMSE)
#' @param doMC A boolean for Parallel computation, default FALSE.
#' @param ncore number of CPU cores for parallel computation, default 1.
#' @param HF if available, a vector containing the optimal bandwidth parameters for each
#' covariate, default NULL.
#' @param H0 A bandwidth value for the starting GWR model, default NULL.
#' @param H2 A bandwidth temporal value for the starting GWR model, default NULL.
#' @param Model Type of Model.
#' @param model A previous model estimated using multiscale_gwr function, default NULL
#' @param get_AICg Boolean, should Global AICc be estimated.
#' @param verbose Boolean, verbose mode.
#' @param control a list of extra control arguments, see MGWRSAR help.
#' @return  Return an object of class mgwrsar
multiscale_gwr<-function(formula,data,coords,kernels='bisq',init='GWR',maxiter=20,nstable=6,tolerance=0.000001,doMC=FALSE,ncore=1,HF=NULL,H0=NULL,H2=NULL,Model=NULL,model=NULL,get_AICg=FALSE,verbose=FALSE,control=list(SE=FALSE,adaptive=TRUE,NN=800,isgcv=FALSE,family=gaussian())){
  family=control$family
  controlv<-control
  if(get_AICg) {
  control$get_Rk=TRUE
  control$get_ts=TRUE
  control$get_s=TRUE
  }
  controlv$criterion='AICc'
  controlv$get_Rk=FALSE
  controlv$get_s=FALSE
  controlv$get_ts=TRUE
  if(is.null(control$TP)) control$TP=1:nrow(coords)
  if(family$family!='gaussian') stop('multiscale_gwr works oly with gaussian family')
  if(is.null(Model)) Model = 'GWR'
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
  if(!('indexG' %in% names(control))) {
    while(sum(duplicated(coords))>0) {
      coords<-jitter(coords,0.001)
    }
    G<-prep_d(coords,control$NN,control$TP)
    controlv$indexG<-control$indexG<-G$indexG
    controlv$dists<-control$dists<-G$dists
  } else if(!is.null(model)){
    if(length(model@indexG)==0) {
      G<-prep_d(coords,control$NN,control$TP)
      controlv$indexG<-control$indexG<-G$indexG
      controlv$dists<-control$dists<-G$dists
    } else {
      controlv$indexG<-control$indexG<-model@indexG
      controlv$dists<-control$dists<-model@dists
    }
  }
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
      if(control$adaptive)  {
        if(kernels[1]=='gauss') lower.bound=2 else lower.bound=2*K
        upper.bound=control$NN-2
        tolerance_GWR=1
      } else {
        imin<-which.min(coords[,1]+coords[,2])
        imax<-which.max(coords[,1]+coords[,2])
        max_dist=sqrt((coords[imin,1]-coords[imax,1])^2+(coords[imin,2]-coords[imax,2])^2)
        lower.bound=0
        upper.bound=max_dist
        tolerance_GWR=tolerance
      }
      res=golden_search_bandwidth(formula = formula, H2=H2,data = data,coords=coords, fixed_vars=NULL,kernels=kernels,Model=Model,control=controlv,lower.bound=lower.bound, upper.bound=upper.bound,tolerance = tolerance_GWR)
      H0= res$minimum
      if(get_AICg) {
        modelGWR<-MGWRSAR(formula = formula, H=c(H0,H2),data = data,coords=coords, fixed_vars=NULL,kernels=kernels,Model=Model,control=control)
      } else  modelGWR<-res$model
    }
    control$get_Rk=FALSE
    BETA=modelGWR@Betav
    H=rep(modelGWR@H,K)
    data$eps<-modelGWR@residuals
    if(verbose) cat('H0=',H0,'\n')
    if(get_AICg) {
      St<-modelGWR@Shat
      Rkk<-Rk<-modelGWR@R_k
      AICc<-aicc_f(modelGWR@residuals,sum(diag(St)),n)
      if(verbose) cat('AICc=',AICc,' RMSE=',modelGWR@RMSE,'\n')
    }
  }
  drmse<-delta_rmse<-sqrt(mean(Y^2))
  isgcv=FALSE
  iter=0
  stable<-rep(1,K)
  if(!is.null(HF)) H=HF
  HBETA=list()
  rmse_list<-c()
  while((abs(delta_rmse)>tolerance | any(stable<nstable)) & iter<=maxiter & length(unique(tail(rmse_list)))<2){
    iter=iter+1
    if(verbose) cat('\n')
    if(verbose) cat(' backfitting ')
    for(k in 1:K){
      var = namesX[k]
      if(verbose)  cat(' ', var,' ')
      data$epst<-data$eps+BETA[,k]*data[,var] %>% as.matrix(ncol=1)
      myformula=as.formula(paste0('epst~',var,'-1'))
      if(is.null(HF)){
        if(stable[k]<nstable){
          res=golden_search_bandwidth(formula = myformula,H2=H2, data = data,coords=coords, kernels=kernels,fixed_vars=NULL,Model=Model,control=controlv,lower.bound=lower.bound, upper.bound=upper.bound,tolerance = tolerance_GWR)
          if(get_AICg) {
            res$model<-MGWRSAR(formula = myformula, H=c(res$minimum,H2),data = data,coords=coords, fixed_vars=NULL,kernels=kernels,Model=Model,control=control)
          }
          modellm<-MGWRSAR(formula = myformula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=H0, Model = 'OLS',control=control)
          AIClm=AIC(lm(myformula,data))
          if(res$objective>AIClm){
            res$minimum=n
            res$objective=AIClm
          }
          H[k]<-h<-res$minimum } else h=H[k]} else h=HF[k]
      if(h==n){
        if(stable[k]>=nstable) {
          modellm<-MGWRSAR(formula = myformula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=H0, Model = 'OLS',control=control)
          AIClm=AIC(lm(myformula,data))
        }
        modelGWR=modellm
        modelGWR@H=c(n,H2)
        modelGWR@Betav=as.matrix(rep(modelGWR@Betac,n),ncol=1)
        modelGWR@residuals=residuals(modellm)
      } else {
        if(stable[k]<nstable & is.null(HF)) modelGWR<-res$model else modelGWR<-MGWRSAR(formula =myformula, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(h,H2), Model = 'multiscale_gwr',control=control)
      }
      BETA[,k]<-modelGWR@Betav
      data$eps=modelGWR@residuals
      if(get_AICg) {
        Sk<-modelGWR@Shat
        Rkk[[k]]<-eigenMapMatMult(Sk,Rk[[k]])  + Sk- eigenMapMatMult(Sk,St)
        St= St + Rkk[[k]] - Rk[[k]]
        new_ts=sum(diag(St))
        AICg<-aicc_f(data$eps,new_ts,n)
        Rk[[k]]<- Rkk[[k]]
      }
      if(iter>1) if(oldH[k]==H[k]) stable[k]=stable[k]+1 else stable[k]=1
    }
    HBETA[[iter]]<-BETA
    rmse_list<-c(rmse_list,rmse)
    oldH<-H
    delta_rmse=(drmse - sqrt(mean(data$eps^2)))/drmse
    drmse= sqrt(mean(data$eps^2))
    if(verbose & get_AICg) cat('\n iter ',iter,' H ',H,' RMSE=',sqrt(mean(data$eps^2)),' delta_rmse ',delta_rmse,' AICg=',AICg,' stable',stable,'\n') else  if(verbose) cat('\n iter ',iter,' H ',H,' RMSE=',sqrt(mean(data$eps^2)),' delta_rmse ',delta_rmse,' stable',stable,'\n')
  }
  if(verbose) cat('Time = ',(proc.time()- start)[3],'\n')
  modelGWR@Model='multiscale_gwr'
  modelGWR@Betav=BETA
  if(get_AICg) modelGWR@AICc=AICg
  modelGWR@residuals=data$eps
  modelGWR@fit=Y-data$eps
  modelGWR@RMSE=sqrt(mean(modelGWR@residuals^2))
  modelGWR@RMSEtp=sqrt(mean(modelGWR@residuals^2))
  modelGWR@H=c(H,H2)
  modelGWR@X=X
  modelGWR@ctime=(proc.time()- start)[3]
  modelGWR
}
