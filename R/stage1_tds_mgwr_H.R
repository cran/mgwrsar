#' stage1_tds_mgwr_H
#' to be documented
#' @usage stage1_tds_mgwr_H(env = parent.frame())
#' @param model  to be documented
#' @noRd
#' @return to be documented
stage1_tds_mgwr_H <- function(env = parent.frame()){
  with(env,{
    i=1
    rmse=sqrt(mean(data$e0^2))
    #browser()
    #data$e0<-e0_init
    #BETA<-BETA_init
    #fit=rowSums(BETA*X)
    #rmse=sqrt(mean(e0_init^2))
    delta_rmse=1
    if(TRUE) while(abs(delta_rmse)>=tol){#& any(stable<2) &
      last_rmseG=rmse
      if(verbose) cat('\n ',i)
      gc()
      myformula_b=as.formula(paste0('e0~-1+',paste0(varying,collapse = '+')))
      controlv<-control
      for(k in varying){
        last_rmse=rmse
        data$e0k<-data$e0+BETA[,k]*X[,k]
        myformula_bk=as.formula(paste0('e0k~-1+',k))
        v=H[k]
        if(v!=max_dist | (v==max_dist & !is.null(H2))) {
          ## metre à jour NN et indexG et dist
          if(control$adaptive){
            if(kernels[1]=='gauss') NNN=round(v*1.5) else NNN=v+2
            controlv$NN=min(NNN,control$NN)
            controlv$indexG=control$indexG[,1:controlv$NN]
            controlv$dists[['dist_s']]= control$dists[['dist_s']][,1:controlv$NN]
          }
          model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(v), Model = 'GWR',control=controlv)
          betav<-model_k@Betav
        } else {
          model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(v), Model = 'OLS',control=controlv)
          betav<-rep(model_k@Betac,n)
        }
        e0=model_k@residuals
        isol=is.na(betav)
        e0[isol] <- data$e0[isol] ## on ne corrige pas les points trop isolé
        BETA[!isol,k]=betav[!isol]
        fit=rowSums(BETA*X)
        data$e0<-Y-fit
      }
      fit=rowSums(BETA*X)
      data$e0<-Y-fit
      rmse=sqrt(mean(data$e0^2))
      delta_rmse=(last_rmse-rmse)/last_rmse
      #if(delta_rmse==0) delta_rmse_0=delta_rmse_0+1
      HBETA[[i+1]]<-BETA
      if(!is.null(TRUEBETA) ){
        for(k in 1:K) HRMSE[i+1,k]=sqrt(mean((TRUEBETA[,k]-BETA[,k])^2))
        HRMSE[i+1,K+1]<-mean(HRMSE[i+1,1:K])
        HRMSE[i+1,K+2]<-opt[k]
        if(get_AIC) HRMSE[i+1,K+3]<-AICc
        HRMSE[i+1,K+4]<-sqrt(mean(data$e0^2))
      }
      #browser()
      if(verbose )  cat('\n delta_rmse ',delta_rmse)

      i=i+1
    }

    ### RETURN MODEL
    modelGWR<-new('mgwrsar')
    modelGWR@mycall<-mycall
    modelGWR@data<-data
    modelGWR@coords<-coords
    modelGWR@X<-modelGWR@XV<-X

    if(verbose)  cat('\n Time = ',(proc.time()- start)[3],'\n')
    modelGWR@formula=formula
    modelGWR@Model=Model
    modelGWR@H=H
    if(control$Type=='GDT') {
      modelGWR@H2=H2
      modelGWR@Z=control$Z
    }
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
    model_stage1<-returned_model<-modelGWR
  })
}
