#' stage1_tds_mgwr
#' to be documented
#' @usage stage1_tds_mgwr(env = parent.frame())
#' @param model  to be documented
#' @noRd
#' @return to be documented
stage1_tds_mgwr <- function(env = parent.frame()){
  with(env,{
    OPT=TRUE
    ### INIT ALGO
    i=1
    if(init_model=='GWR') {
      i=max(which.min(abs(V-h_init))-1,1)
      ts<-model_lm0@tS

    }
    else {
      new_ts<-ts<-model_lm0@tS
      new_TS<-model_lm0@TS
    }

    if(verbose) {
      if(!is.null(model_stage1)) {
        cat('Starting from a previous model : \n')#
        summary(model_stage1)
      } else  cat('\n i=',i,' Starting model \n')
    }
    lvarying=length(varying)
    rmse<-sqrt(mean(data$e0^2))
    if(verbose) cat('\n\n########################################################################\n STAGE 1 : find unique bandwidth for each covariate  \n using Top Down Scale Approach with backfitting ...\n########################################################################\n')
    stable<-rep(0,length(varying))
    up<-rep(V[max(i-1,1)],lvarying)
    last_opt<-opt<-rep(V[i],lvarying)
    ddown<-down<-rep(V[min(i+1,length(V))],lvarying)
    names(stable)<-names(up)<- names(opt)<- names(down)<- names(ddown)<-varying
    if(get_AIC) {
      BestCrit<-lastCrit <-lastCrit2<-10^6
      HTS[1]<-sum(TSik)
      AIC_deep<-NA
      prev_AICc<-AICc+10^6
    }
    if(doMC) {
      registerDoParallel(cores=ncore)
    } else registerDoSEQ()
    bestBETA=BETA
    delta_rmse=1
    spacestable=TRUE
    while(abs(delta_rmse)>=tol){
      last_rmseG=rmse
      if(verbose) cat('\n ',i)
      if(get_AIC) {
        prev_AICc=AICc
        BestCrit=lastCrit
        St=S
      }
      varyingT<-varying
      if(i==browser) browser()
      gc()
      myformula_b=as.formula(paste0('e0~-1+',paste0(varyingT,collapse = '+')))
      controlv<-control
      last_opt=opt
      for(k in varying){
        if(verbose) cat(' ',k)
        last_rmse=rmse
        if(opt[k]<max_dist) up[k]<-tail(V[V>opt[k]],1) else up[k]<-max_dist
        if(opt[k]>min_dist) {
          if(sum(V<opt[k])>0) down[k]<-head(V[V<opt[k]],1) else down[k]<-min_dist
          if(sum(V<min(down))>0) {
            ddown[k]<-head(V[V<min(down)],1)
          } else  ddown[k]<-down[k]
        } else ddown[k]<-down[k]<-min_dist
        data$e0k<-data$e0+BETA[,k]*X[,k]
        myformula_bk=as.formula(paste0('e0k~-1+',k))
        vks=c(up[k],opt[k],down[k],ddown[k])
        if(any(stable<2)){
          res<-foreach(v =unique(vks),.combine="rbind",.inorder=FALSE)  %do% {
            if(v!=max_dist | (v==max_dist & !is.null(H2))) {
              if(control$adaptive[1]){
                if(kernels[1]=='gauss') NNN=round(v*1.5) else NNN=v+2
                controlv$NN=min(NNN,control$NN)
                controlv$indexG=control$indexG[,1:controlv$NN]
                controlv$dists[['dist_s']]= control$dists[['dist_s']][,1:controlv$NN]
              }
              model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(v, NULL), Model = 'GWR',control=controlv)
              Betav<-model_k@Betav
            } else {
              model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(v,NULL), Model = 'OLS',control=controlv)
              Betav<-rep(model_k@Betac,n)
              model_k@AICc=model_k@AIC
            }
            list(AICc=model_k@AICc,betav=Betav,e0=residuals(model_k),vk=v,TS=as.numeric(model_k@TS),S=model_k@Shat)
          }
          mybest= which.min(res[,'AICc'])
          opt[k]<-unlist(res[mybest,'vk'])
          if(any(stable>0) & !spacestable) {
            spacestable=TRUE
          }
          e0=unlist(res[mybest,'e0'])
          betav<-unlist(res[mybest,'betav'])
          isol=is.na(betav)
          e0[isol] <- data$e0[isol]
          if(get_AIC) {
            TSik[!isol,k]<-unlist(res[mybest,'TS'])[!isol]
            Sk<-unlist(res[mybest,'S'][[1]])
          }
        } else {
          if(opt[k]<n-2)  model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=c(opt[k],NULL), Model = 'GWR',control=controlv) else model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=c(opt[k],NULL), Model = 'OLS',control=controlv)
          e0=residuals(model_k)
          betav<-model_k@Betav
          isol=is.na(betav)
          e0[isol] <- data$e0[isol]
          if(get_AIC ) {
            TSik[!isol,k]<-model_k@TS
            Sk<-model_k@Shat
          }
        }
        BETA[!isol,k]=betav[!isol]
          if(get_AIC) {
            Rkk[[k]]<-eigenMapMatMult(Sk,Rk[[k]])  + Sk- eigenMapMatMult(Sk,St)
            St= St + Rkk[[k]] - Rk[[k]]
            new_TS=diag(St)
            new_ts=sum(new_TS)
            if(verbose & get_AIC)  cat('\n AICc : ',  aicc_f(e0[idx_init],new_ts,n_time))
          }
          fit=rowSums(BETA*X)
          data$e0<-Y-fit
        rmse=sqrt(mean(data$e0^2))
      }
      for(k in fixed_vars){
        last_rmse=rmse
        data$e0k<-data$e0+BETA[,k]*X[,k]
        myformula_bk=as.formula(paste0('e0k~-1+',k))
        model_tds<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=NULL, Model = 'OLS',control=controlv)
          BETA[,k]=model_tds@Betac
          if(get_AIC) {
            TSik[,k]<-unlist(res[mybest,'TS'])
            Sk<-unlist(res[mybest,'S'][[1]])
            Rkk[[k]]<-eigenMapMatMult(Sk,Rk[[k]])  + Sk- eigenMapMatMult(Sk,St)
            St= St + Rkk[[k]] - Rk[[k]]
            new_TS=diag(St)
            new_ts=sum(new_TS)
          }
          fit=rowSums(BETA*X)
          data$e0<-Y-fit
        rmse=sqrt(mean(data$e0^2))
      }
      if(get_AIC) {
        AICc<-aicc_f(data$e0[idx_init],new_ts,n_time)
      }
      bestBETA=BETA
      if(get_AIC) {
        S=St
        for(k in varying) Rk[[k]]<-Rkk[[k]]
      }
      stable[opt!=last_opt]<-0
      stable[opt==last_opt]<-stable[opt==last_opt]+1
      fit=rowSums(BETA*X)
      data$e0<-Y-fit
      rmse=sqrt(mean(data$e0^2))
      delta_rmse=(last_rmseG-rmse)/last_rmseG
      if(get_AIC) {
        HTS[[i+1]]<-new_TS
        HAICc<-c(HAICc,AICc)
      }
      HBETA[[i+1]]<-BETA
      if(!is.null(TRUEBETA) ){
        for(k in 1:K) HRMSE[i+1,k]=sqrt(mean((TRUEBETA[,k]-BETA[,k])^2))
        HRMSE[i+1,K+1]<-mean(HRMSE[i+1,1:K])
        HRMSE[i+1,K+2]<-opt[k]
        if(get_AIC) HRMSE[i+1,K+3]<-AICc
        HRMSE[i+1,K+4]<-sqrt(mean(data$e0^2))
      }
      if(verbose & get_AIC)  cat('\n AICc : ', AICc, ' stable ',stable,' delta_rmse ',delta_rmse)
      if(verbose & !get_AIC)  cat('\n stable ',stable,' delta_rmse ',delta_rmse)
      if(verbose) cat('\n\n H=',opt)
      i=i+1
    }
    H<-opt
    names(H)<-varying
    mybest<-i
    if(get_AIC) {
      AICc<-HAICc[mybest]
      TS=as.numeric(unlist(HTS[mybest]))
      tS=sum(TS)
    }
    BETA<-HBETA[[mybest]]
    if(!is.null(TRUEBETA) ){ HRMSE<-HRMSE[1:mybest,]}
    HBETA<-HBETA[1:mybest]
    fit=rowSums(BETA*X)
    data$e0<-Y-fit
    ### RETURN MODEL
    modelGWR<-new('mgwrsar')
    modelGWR@mycall<-mycall
    modelGWR@data<-data
    modelGWR@coords<-as.matrix(coords)
    modelGWR@X<-modelGWR@XV<-X
    if(verbose)  cat('\n Time = ',(proc.time()- start)[3],'\n')
    modelGWR@formula=formula
    modelGWR@Model=Model
    modelGWR@H=H
    modelGWR@fixed_vars<-unique(c(as.character(fixed_vars),names(modelGWR@H)[which(modelGWR@H==max_dist)]))
    modelGWR@Betav=BETA
    modelGWR@fit=rowSums(BETA*X)
    modelGWR@residuals=Y-modelGWR@fit
    modelGWR@RMSE=sqrt(mean(modelGWR@residuals^2))
    modelGWR@Type=control$Type
    modelGWR@kernels=kernels
    modelGWR@adaptive=control$adaptive
    modelGWR@TP=1:n
    modelGWR@doMC=control_tds$doMC
    modelGWR@ncore=control_tds$ncore
    modelGWR@NN=control$NN
    modelGWR@alpha=control$alpha
    modelGWR@V=V
    modelGWR@Y=Y
    modelGWR@ctime <- (proc.time()- start)[3]
    modelGWR@HBETA<-HBETA
    if(get_AIC) {
      if(!is.null(S)) modelGWR@Shat=S
      modelGWR@tS=tS
      modelGWR@TS=TS
      modelGWR@edf <- n-tS
      modelGWR@AICc <-AICc
      for(k in varying){
        modelGWR@R_k[[k]]<- Rk[[k]]
      }
    }
    if(!is.null(TRUEBETA)) {
      modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[,1]),]
    }
    modelGWR@G=G
    control_tds$model_stage1<-model_stage1<-returned_model<-modelGWR
  })
}
