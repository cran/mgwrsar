#' atds_gwr Top-Down Scaling approach of GWR
#'
#' This function performs a Geographically Weighted Regression (GWR) using
#' a top-down scaling approach, adjusting GWR coefficients with a progressively
#' decreasing bandwidth as long as the AICc criterion improves.
#'
#' @usage atds_gwr(formula,data,coords,kernels='triangle',fixed_vars=NULL,
#' control_tds=list(nns=30),control=list(adaptive=TRUE,verbose=FALSE))
#' @param formula  a formula.
#' @param data a dataframe.
#' @param coords default NULL, a dataframe or a matrix with coordinates
#' @param kernels A vector containing the kernel types. Possible types:
#' triangle ("triangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane").
#' @param fixed_vars a vector with the names of spatiallay constant
#' coefficient for mixed model. All other variables present in formula
#' are supposed to be spatially varying. If empty or NULL (default),
#' all variables in formula are supposed to be spatially varying.
#' @param control list of extra control arguments for MGWRSAR wrapper - see MGWRSAR Help
#' @param control_tds list of extra control arguments for tds_mgwr model - see tds_gwr Help
#' @seealso  tds_mgwr, gwr_multiscale, MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar.
atds_gwr<-function(formula,data,coords,kernels='triangle',fixed_vars=NULL,control_tds=list(nns=30),control=list(adaptive=TRUE,verbose=FALSE)){
  Model='atds_gwr'
  #criteria='autre'
  start<-proc.time()
  init_param_tds() #done ?
  if(is.null(control_tds$V)) built_Vseq() #done ?
  if(!is.null(TRUEBETA)) init_RMSE_history()   #done ?

  ### INIT ALGO
  e0=residuals(model_lm0)
  ds0<-TS<-model_lm0@TS
  ds0=diag(S)

  i=1
  ts=1
  #ndown=1
  tocorrect<-1:n
  BestCrit<-lastCrit <-lastCrit2<-10^6
  if(verbose) cat('\n i=',i,' Starting model \n')#
  CONTINUE=TRUE
  lvarying=length(varying)
  if(length(varying)==0) {
    returned_model<- model_lm0

  } else {
 # if(browser==0) browser()

    ### LOOP ALGO
    while( (CONTINUE & i<length(V))){ # | (any(H<=V[i] & control$isgcv))
      #browser()
      #if(all(H>V[i])) break
      varyingT<-varying
      if(i==browser) browser()
      gc()
      if(length(varying>1) & !is.null(H)) for(k in varying){
          if(H[k]>V[i]) {
            varyingT<-setdiff(varyingT,k)
          }
      }
      vk=V[i]

      if(i>1) vkm1=V[i-1] else {
        if(control$adaptive[1])  vkm1=control$NN-2 else {
          imin<-which.min(coords[,1]+coords[,2])
          imax<-which.max(coords[,1]+coords[,2])
          max_dist=sqrt((coords[imin,1]-coords[imax,1])^2+(coords[imin,2]-coords[imax,2])^2)
          vkm1=max_dist
        }

      }
      myformula_b=as.formula(paste0('e0~-1+',paste0(varyingT,collapse = '+')))
      controlv<-control
      if(control$adaptive[1]) controlv$NN<-vk+2 else controlv$NN=nrow(coords)
      controlv$get_s=TRUE
      if(verbose) cat('\n GWR Correction of varying coefficients with v=',vk, ' for varyings coefficients: ', paste(varyingT,collapse=' '))
      modelGWR<-MGWRSAR(formula = myformula_b, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(vk,control_tds$H2),Model = 'GWR',control=controlv)
      ##### diagnostic
      e1=residuals(modelGWR)

      #browser()


        S1<-S+modelGWR@Shat-eigenMapMatMult(modelGWR@Shat,S)
        ds1=diag(S1)
        df_true=sum(diag(ds1))
        tocorrect<-1:n
        S[tocorrect,]=S1[tocorrect,]
        ds0[tocorrect]=ds1[tocorrect]
      BETA[tocorrect,colnames(modelGWR@Betav)]=BETA[tocorrect,colnames(modelGWR@Betav)]+modelGWR@Betav[tocorrect,]
      fit=rowSums(BETA*X)
      e0<-data$e0<-Y-fit

      if(length(fixed_vars)>0) {
        if(verbose) cat('\n LM Correction of non varying coefficients: ', paste(fixed_vars,collapse=' '))
        formula_constant=as.formula(paste0('e0~-1+',paste0(fixed_vars,collapse = '+')))
        model_cor_constant=lm(formula_constant,data=data)
        for(k in fixed_vars) {
          BETA[,k]=BETA[,k]+coef(model_cor_constant)[k]
        }
        fit=rowSums(BETA*X)
        data$e0<-Y-fit
      }
        global_ts<-df_true
        local_ts<-sum(ds0)
      lastCrit<-AICc<-n*log(mean(data$e0^2))+n*log(2*pi)+n*(n+global_ts)/(n-2-global_ts)
      if(verbose) cat('\n i=',i,' v=',vk,' AICc = ',AICc,' df_true = ',df_true,'\n')

      if((BestCrit-lastCrit)/abs(BestCrit) >=tol ) {
        if(sum(tocorrect)<=round(n*0.333) & i>1) CONTINUE=FALSE else {
        mybestBETA=BETA
        ## best diagnostic
        S_best=S
        TS_best=diag(S)
        tS_best=sum(TS_best)
        AIC_opt<-AICc
        HOPT[varyingT]<-vk
        BestCrit=lastCrit
        }

      }
      if(length(varying)==0) CONTINUE=FALSE
      HBETA[[i+1]]<-BETA
      if(!is.null(TRUEBETA) ){
        for(k in 1:K) HRMSE[i+1,k]=sqrt(mean((TRUEBETA[,k]-BETA[,k])^2))
        HRMSE[i+1,K+1]<-mean(HRMSE[i+1,1:K])
        HRMSE[i+1,K+2]<-vk
        HRMSE[i+1,K+3]<-AICc
        HRMSE[i+1,K+4]<-sum(data$e0^2)
        HRMSE[i+1,K+5]<-global_ts
      }
      i=i+1
    }
    ### RETURN MODEL
    if(verbose)  cat(' Time = ',(proc.time()- start)[3],'\n')
    modelGWR@Model='atds_gwr'
    modelGWR@fixed_vars<-as.character(fixed_vars)
    if(control$isgcv) {
      modelGWR@Betav=BETA
      modelGWR@H=H
      modelGWR@AICc <- AICc
      } else  {
        modelGWR@H=HOPT
        modelGWR@Betav=mybestBETA
        modelGWR@AICc <- AIC_opt
      }
    modelGWR@fit=rowSums(modelGWR@Betav*X)
    modelGWR@residuals=Y-modelGWR@fit
    modelGWR@RMSE=sqrt(mean(modelGWR@residuals^2))
    modelGWR@V=c(n,V)
    modelGWR@X=X
    modelGWR@Y=Y
    modelGWR@ctime <- (proc.time()- start)[3]
    modelGWR@HBETA<-HBETA
    if(!is.null(TRUEBETA)) {
      colnames(HRMSE)<-c(paste0('RMSE_',namesX),'meanRMSE','v','AICc','SSR','TS')
      modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[,1]),]
    }
    modelGWR@G<-list(indexG= control$indexG,dists=control$dists)
    modelGWR@Shat<-S_best
    modelGWR@TS<-TS_best
    modelGWR@tS<-tS_best
    modelGWR@edf <- n-tS_best
    returned_model<-modelGWR
  }
  returned_model
}

