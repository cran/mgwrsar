#' stage2_atds_mgwr
#' to be documented
#' @usage stage2_atds_mgwr(env = parent.frame())
#' @param model  to be documented
#' @noRd
#' @return to be documented
stage2_atds_mgwr <- function(env = parent.frame()){
  with(env,{
    if(browser==-2) browser()
    if(verbose) cat('Start stage 2 : \n')
    nrounds=control_tds$nrounds
    if(!is.null(control_tds$model_stage1)) model_stage1<-control_tds$model_stage1
    if(length(model_stage1@R_k)==0) stop('Starting model_stage1 must be runned with get_AIC=TRUE')
    if(verbose) summary(model_stage1)
    control_tds$V<-model_stage1@V
    varying<-setdiff(namesX,fixed_vars)
    controlv<-control
    controlv$get_s=TRUE
    if(!is.null(control_tds$model_stage1)){
      HBETA=model_stage1@HBETA
      HRMSE=model_stage1@HRMSE
      i<-nrow(HRMSE)
      HRMSE<-rbind(HRMSE,matrix(NA,ncol=ncol(HRMSE),nrow=nrounds+1))
      BETA<-model_stage1@Betav
      AICc<-model_stage1@AICc
      varying<-varying[order(model_stage1@H,decreasing = T)]
      fixed_vars<-model_stage1@fixed_vars
      varyingT<-setdiff(varying,fixed_vars)
      data$e0<-residuals(model_stage1)
      TS=model_stage1@TS
      H=model_stage1@H
      G<-model_stage1@G
      control$indexG=G$indexG
      control$dists=G$dists
    }
    if(verbose) cat('\n\n########################################################################\n STAGE 2 :  use multiple/adaptive bandwidth for each covariate \n using Top Down Scale Approach with backfitting ...\n########################################################################\n')
    HH<-list()
    nround=1
    n<-nrow(model_stage1@Betav)
    if(control_tds$get_AIC){
      Rkk<-Rk<-model_stage1@R_k
      S<-model_stage1@Shat
      TS<-model_stage1@TS
      tS=model_stage1@tS
      last_AICck<-last_AICc<-AICc+1
    } else {
      AICc=0;
      last_AICck<-last_AICc<-1
      tS<-ts<-S<-NULL
    }
    rmse=sqrt(mean(data$e0^2))+10^6
    n_updated=1
    while(nround<=nrounds & n_updated>0){
      last_AICck<-last_AICc<-AICc
      mybestbeta<-BETA
      mybestAICc<-AICc
      mybesttS<-tS
      mybestTS<-TS
      mybestShat=S
      mybestedf<-n-tS
      n_updated=0
      for(k in namesX){
        rmse<-sqrt(mean(data$e0^2))
        control_tds2=control_tds
        control_tds2$H<-0
        names(control_tds2$H)=k
        control_tds2$V<- control_tds2$V[-1]
        control_tds2$model_stage1=NULL
        control_tds2$TRUEBETA=TRUEBETA
        if(!is.null(TRUEBETA)) {
          colnames(control_tds2$TRUEBETA)<-namesX
          control_tds2$TRUEBETA<-matrix(control_tds2$TRUEBETA[,k],ncol=1)
        }
        data$e0k<-data$e0+BETA[,k]*X[,k]
        myformula_bk=as.formula(paste0('e0k~-1+',k))
        if(k %in% fixed_vars) {
          model_tds<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=NULL, Model = 'OLS',control=controlv)
          HH[[k]]<-n
          model_tds@Betav<-as.matrix(rep(model_tds@Betac,n),ncol=1)
        } else {
          model_tds<-atds_gwr(formula=myformula_bk,data=data,coords=coords,kernels=kernels,fixed_vars=NULL,control_tds=control_tds2,control=controlv)
        }
        betav<-model_tds@Betav
        e0=residuals(model_tds)
        idxt=!is.na(betav)
        if(control_tds$get_AIC){
          Sk<-model_tds@Shat
          Rkk[[k]]<-eigenMapMatMult(Sk,Rk[[k]]) + Sk- eigenMapMatMult(Sk,S)
          St= S + Rkk[[k]] - Rk[[k]]
          new_ts=sum(diag(St))
          model_tds@AICc<-n*log(sum(e0[idxt]^2)/n_time)+n_time*log(2*pi)+n_time*(n_time+new_ts)/(n_time-2-new_ts)
        }
        if(nround<nrounds) fullupdate=FALSE else fullupdate=TRUE
        if( model_tds@AICc<last_AICck |  sqrt(mean(e0^2))<rmse ) {
          if(verbose) cat(paste0(' ',k,' updated; '))
          n_updated=n_updated+1
          BETA[idxt,k]=betav[idxt]
          last_AICck<-AICc<-model_tds@AICc
          if(!(k %in% fixed_vars)) HH[[k]]<-model_tds@V[model_tds@V>=model_tds@H]
          fit=rowSums(BETA*X)
          data$e0<-Y-fit
          if(control_tds$get_AIC){
            S=St
            TS=diag(S)
            tS=sum(TS)
            Rk[[k]]<-Rkk[[k]]
          }
        } else {
          if(verbose)  cat(paste0(' ',k,' not updated; '))
        }
      }
      if(!(last_AICc-AICc)/abs(AICc)>tol & nrounds>=10){
        BETA<-mybestbeta
        AICc<-mybestAICc
        tS<-mybesttS
        TS<-mybestTS
        S<-mybestShat
      }
      HBETA[[length(HBETA)+1]]<-BETA
      if(!is.null(TRUEBETA) ){
        for(k in 1:K) HRMSE[i+1,k]=sqrt(mean((TRUEBETA[,k]-BETA[,k])^2))
        HRMSE[i+1,K+1]<-mean(HRMSE[i+1,1:K])
        cat('\n BETA RMSE = ',mean(HRMSE[i+1,1:K]),'\n')
        HRMSE[i+1,K+2]<-min(unlist(HH))
        HRMSE[i+1,K+3]<-AICc
        HRMSE[i+1,K+4]<-sqrt(mean(data$e0^2))
      }
      nround= nround+1
      i=i+1
      if( !is.null(TRUEBETA) & verbose) cat('\n nround ',nround-1,' mean beta rmse ',HRMSE[i+1,K+1],'  rmse ',sqrt(mean(data$e0^2)),'\n')
      if( is.null(TRUEBETA)  & verbose ) cat('\n nround ',nround-1,  ' rmse ',sqrt(mean(data$e0^2)),'\n')
    }
    if(verbose) for(k in varying){
      cat('\n',k,' : ',unlist(HH[[k]] ),'\n')
    }
    ### RETURN MODEL
    modelGWR<-model_stage1
    modelGWR@mycall<-mycall
    modelGWR@data<-data
    modelGWR@coords<-as.matrix(coords)
    modelGWR@X<-modelGWR@XV<-X
    modelGWR@formula=formula
    if(verbose)  cat('\n Time = ',(proc.time()- start)[3],'\n')
    modelGWR@Model=Model
    modelGWR@fixed_vars<-as.character(fixed_vars)
    modelGWR@Betav=BETA
    modelGWR@fit=rowSums(BETA*X)
    modelGWR@residuals=Y-modelGWR@fit
    modelGWR@RMSE=sqrt(mean(modelGWR@residuals^2))
    myH<-H
    for(k in varying){
      if(!is.null(HH[[k]])) {
        myH[k]<-min(unlist(HH[[k]]))
      }
    }
    modelGWR@H=myH
    if(model_stage1@Type=='GDT'){
      modelGWR@Type=model_stage1@Type
      modelGWR@kernels=model_stage1@kernels
      modelGWR@H2=model_stage1@H2
      modelGWR@adaptive=model_stage1@adaptive
    } else {
      modelGWR@Type=control$Type
      modelGWR@kernels=kernels
      modelGWR@adaptive=control$adaptive
    }
    modelGWR@TP=1:n
    modelGWR@doMC=control_tds$doMC
    modelGWR@ncore=control_tds$ncore
    modelGWR@TP=1:n
    modelGWR@V=c(n,V)
    modelGWR@Y=Y
    modelGWR@ctime <-  (proc.time()- start)[3]
    modelGWR@HBETA<-HBETA
    if(!is.null(TRUEBETA)) {
      modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[,1]),]
    }
    if(control_tds$get_AIC) {
      modelGWR@AICc <-AICc
      modelGWR@tS=tS
      modelGWR@TS=TS
      modelGWR@Shat=S
      modelGWR@edf<-n-tS
    }
    returned_model<-modelGWR
  })
}


