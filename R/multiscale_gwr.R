#' multiscale_gwr
#' to be documented (experimental)
#' @usage multiscale_gwr(formula, data,coord, fixed_vars=NULL,Model=c('GWR'),
#' kernels=c('bisq'),maxiter=10,control=list(),fd=c(0.6,0.9,0.99,1,1.01,1.1,1.4),
#' H=NULL,BETA=NULL,grad=FALSE,nstable=3,nfull=2,init='GWR')
#' @param formula  a formula.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coord default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param Model character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' explanation.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"), gaussian
#' ("gauss")) .
#' @param H vector containing the bandwidth parameters for each covariate.
#' @param maxiter maximum number of backfiting iteration
#' @param BETA to be documented (experimental)
#' @param nstable to be documented (experimental)
#' @param init to be documented (experimental)
#' @param fd to be documented (experimental)
#' @param grad to be documented (experimental)
#' @param nfull to be documented (experimental)
#' @param control list of extra control arguments for MGWRSAR wrapper - see Details below

multiscale_gwr<-function(formula, data,coord, fixed_vars=NULL,Model=c('GWR'),kernels=c('bisq'),maxiter=10,control=list(),fd=c(0.6,0.9,0.99,1,1.01,1.1,1.4),H=NULL,BETA=NULL,grad=FALSE,nstable=3,nfull=2,init='GWR'){
  # Stage 1 global bandwidth
  #control$verbose=0
  controlGCV=control
  controlGCV$isgcv=TRUE
  RES<-list()
  find_end_cycle<-function(v,prof){
    cycling=FALSE
    myprof=NA
    fin=length(v)
    for(l in 2:prof){
      mycycle=tail(v,l)
      for(j in 0:(l-1)){
        if(all(tail(v[-((fin-j):fin) ],l)==mycycle)) {
          cycling=TRUE
          myprof=l
        }
      }
    }
    list(iscycling=cycling,prof=myprof)
  }
  mytab<-bandwidths_mgwrsar(formula = formula, data = data,coord=coord, fixed_vars=fixed_vars,Models=Model,candidates_Kernels=kernels,control=control)
  initmodelGCV<-MGWRSAR(formula = formula, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=mytab[[1]]$model$H, Model = Model,control=controlGCV)
  initmodel<-MGWRSAR(formula = formula, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=mytab[[1]]$model$H, Model = Model,control=control)
  RES[[1]]<-list(Betav=initmodelGCV$Betav,residuals=initmodelGCV$residuals,fit=initmodelGCV$fit,RMSE=initmodelGCV$RMSE,H=initmodelGCV$H,CV=initmodelGCV$RMSE)

  if(!is.null(BETA)) cat('\n STARTING RMSE_beta=',apply(mytab[[1]]$model$Betav-BETA,2,function(x) sqrt(mean(x^2))),'\n')
  ##### Stage 2 backfit
  # INIT
  namesX<-attr(terms(formula),"term.labels")
  if(attr(terms(formula),'intercept')==1) {
    intercept=TRUE
    namesX=c('Intercept',namesX)
    data$Intercept<-rep(1,nrow(data))
  }
  m=length(namesX)
  for(i in 1:(m)){
    RES[[i]]<-RES[[1]]
  }
  started<-convcrit<-FALSE

  iter=2
  stable=rep(0,m)
  RMSEs<-NA
  CVs<-rmse<-RES[[1]]$RMSE
  H_p<-matrix(NA,ncol=maxiter+1,nrow=m)
  H_p[,1]<-RES[[1]]$H

  ## FIND H using LOOCV
  cat('\n##############################\n')
  cat('########## FIND H using LOOCV \n')
  cat('##############################\n')
  if(is.null(H)){
    instable=1:m
    while(length(instable)>0){
      instable=which(stable<nstable)
      for(i in 1:m){
        control$isgcv=FALSE
        myformula_b<-as.formula(paste0('Ybackfit~-1+',namesX[i]))
        cat('\n VARIABLE : ',namesX[i] ,' stable since ',stable[i],' iterations')
        if(i==1) j=m else j=i-1
        data$Ybackfit<-RES[[j]]$residuals+RES[[i]]$Betav[,namesX[i]]*data[,namesX[i]]
        ## soit H connu soit non
        if(i %in% instable){
          if(grad & iter>nfull) {
            RES[[i]]<-bandwidth_backfit(formula = myformula_b, data = data,coord=coord, H= RES[[i]]$H,fixed_vars=fixed_vars,Model=Model,kernels=kernels,control=control,fd=fd,isgcv=TRUE)
            H_p[i,iter]<-RES[[i]]$H
          } else {
            myres<-bandwidths_mgwrsar(formula = myformula_b, data = data,coord=coord, fixed_vars=fixed_vars,Models=Model,candidates_Kernels=kernels,control=control,control_search=list())
            model<-MGWRSAR(formula = myformula_b, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=myres[[1]]$model$H, Model = Model,control=controlGCV)
            RES[[i]]<-list(Betav=model$Betav,fit=model$fit,residuals=model$residuals,CV=model$RMSE,RMSE=model$RMSE,H=model$H)
            H_p[i,iter]<-RES[[i]]$H
          }
          cat('\n new value ',RES[[i]]$H)
          if(H_p[i,iter]==H_p[i,iter-1]) stable[i]=stable[i]+1
          if(iter>5 & stable<nstable) {
            test_cylcing<-find_end_cycle(H_p[i,1:iter],4)
            if(test_cylcing$iscycling) {
              H_p[i,iter]=ceiling(mean(H_p[i,((iter-prof+1):iter)]))
              stable[i]=nstable
            }
          }
        } else {
          H_p[i,iter]<- H_p[i,iter-1]
          modelGCV<-MGWRSAR(formula = myformula_b, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=H_p[i,iter], Model = Model,control=controlGCV)
          RES[[i]]<-list(Betav=modelGCV$Betav,fit=modelGCV$fit,residuals=modelGCV$residuals,CV=modelGCV$RMSE,RMSE=modelGCV$RMSE,H=H_p[i,iter])
        }
      }
      deltaCVp=(rmse-RES[[m]]$RMSE)/RES[[m]]$RMSE
      rmse=RES[[m]]$RMSE
      CVs<-c(CVs,rmse)
      RMSEs<-c(RMSEs,NA)
      if(any(stable==nstable) & all(H_p[,iter]==H_p[,iter-2])) {
       # browser()
        stable[1:m]=nstable
        for(i in 1:m) H_p[i,iter]=round((H_p[i,iter]+H_p[i,iter-1])/2)
      }
      cat('\n bandwidths : ',H_p[,iter],' deltaCVp =',deltaCVp,'\n ')
      iter=iter+1
    }
  } else   {
    H_p[,iter]<-H
    iter=iter+1
  }
  cat('\n##############################\n')
  cat('########## BETA estimation    \n')
  cat('##############################\n')
  RES<-list()
  if(init=='0'){
    initmodel$Betav=matrix(0,nrow=nrow(data),ncol=ncol(initmodel$Betav),byrow=T)
    colnames(initmodel$Betav)<-namesX
    initmodel$residuals=initmodel$Y #residuals(lmmodel)
    initmodel$fit=rep(0,nrow(data)) #fitted(lmmodel)
    initmodel$RMSE=sqrt(mean(initmodel$residuals^2))
  } else if(init=='lm'){
    lmmodel=lm(formula = formula, data = data)
    initmodel$Betav=matrix(coefficients(lmmodel),nrow=nrow(data),ncol=ncol(initmodel$Betav),byrow=T)
    colnames(initmodel$Betav)<-namesX
    initmodel$residuals=residuals(lmmodel)
    initmodel$fit=fitted(lmmodel)
    initmodel$RMSE=sqrt(mean(initmodel$residuals^2))
  }

  RES[[1]]<-list(Betav=initmodel$Betav,residuals=initmodel$residuals,fit=initmodel$fit,RMSE=initmodel$RMSE,H=initmodel$H,CV=initmodel$RMSE)
  for(i in 1:(m)){
  RES[[i]]<-RES[[1]]
  }
  SRES<-RES
  started<-convcrit<-FALSE
  rmse=initmodel$RMSE
  iterbeta=1
  while(!convcrit & iter<maxiter){
    for(i in 1:m){
      control$isgcv=FALSE
      myformula_b<-as.formula(paste0('Ybackfit~-1+',namesX[i]))
      if(i==1) j=m else j=i-1
      data$Ybackfit<-RES[[j]]$residuals+RES[[i]]$Betav[,namesX[i]]*data[,namesX[i]]
      H_p[i,iter]<- H_p[i,iter-1]
      model<-MGWRSAR(formula = myformula_b, data = data,coord=coord, fixed_vars=fixed_vars,kernels=kernels,H=H_p[i,iter], Model = Model,control=control)
      ### taux apprentissage ??
      #model$Betav=nu*RES[[i]]$Betav[,namesX[i]]+(1-nu)*model$Betav
      model$fit=model$Betav*data[,namesX[i]]
      model$residuals=data$Ybackfit- model$fit
      model$RMSE=sqrt(mean(model$residuals^2))
      RES[[i]]<-list(Betav=model$Betav,fit=model$fit,residuals=model$residuals,CV=model$RMSE,RMSE=model$RMSE,H=H_p[i,iter])
    }
    deltaRMSEp=(rmse-RES[[m]]$RMSE)/RES[[m]]$RMSE
    if(deltaRMSEp>0)  SRES<-RES
    cat('\n deltaRMSEp =',deltaRMSEp,'\n ')
    if(deltaRMSEp<0.00001 & iterbeta>2) {
      convcrit=TRUE;
      cat('\n Delta RMSE <0.00001\n')
    }
    rmse=RES[[m]]$RMSE
    CVs<-c(CVs,NA)
    RMSEs<-c(RMSEs,rmse)
    iter=iter+1
    iterbeta=iterbeta+1
  }
  cat('\n##############################\n')
  cat('########## output model   \n')
  cat('##############################\n')

  mymodel=mytab[[1]]$model

  for(i in 1:m){
    mymodel$Betav[,i]<-SRES[[i]]$Betav[,1]
  }
  mymodel$fit=SRES[[i]]$fit
  mymodel$residuals=SRES[[i]]$residuals
  mymodel$RMSE=sqrt(mean(mymodel$residuals^2))
  mymodel$Hs=H_p[,1:(iter-1)]
  mymodel$H=H_p[,(iter-1)]
  mymodel$iter=iter
  mymodel$RMSE_iter=RMSEs
  mymodel$CV_iter=CVs
  mymodel$Model<-'multiscale_GWR'
  mymodel
}
