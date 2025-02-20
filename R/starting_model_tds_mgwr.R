#' starting_model_tds_mgwr
#' to be documented
#' @usage starting_model_atds_gwr(env = parent.frame())
#' @param env an environment
#' @return an environment.
#' @noRd
starting_model_tds_mgwr<-function(env = parent.frame()){
  with(env,{
    #browser()
    if(init_model=='GWR'){
      controlv<-control
      #controlv$adaptive=T
      if(controlv$adaptive[1]){
        lower.bound=length(namesX)+2
        upper.bound=controlv$NN
      } else {
        # imin<-which.min(coords[,1]+coords[,2])
        # imax<-which.max(coords[,1]+coords[,2])
        # max_dist=sqrt((coords[imin,1]-coords[imax,1])^2+(coords[imin,2]-coords[imax,2])^2)
        lower.bound=0
        upper.bound=max_dist
      }
      res<-golden_search_bandwidth(formula = formula, H2=H2,data = data,coords=coords,fixed_vars=fixed_vars,kernels=kernels, Model = 'GWR',control=controlv,lower.bound=lower.bound,upper.bound=upper.bound)
      model_lm0<-res$model
      h_init=res$minimum
      BETA=model_lm0@Betav
      } else if(init_model=='OLS'){
      model_lm0<-MGWRSAR(formula = formula, data = data,coords=coords,fixed_vars=NULL, Model = 'OLS',H=NULL,kernels=NULL,control=list())
      BETA=matrix(coef(model_lm0)$Betac,byrow=TRUE,nrow=,nrow(data),ncol=K)
      if(control$adaptive[1]) h_init=n else h_init=max_dist
    } else stop('init_model should be OLS or GWR')
    #model_lm0=lm(formula,data=data)
    colnames(BETA)<-namesX
    data$e0<-residuals(model_lm0)
    idx_init<-idx<-1:n
    drop_init=c()
    n_time<-n
    # S or TS init

    if(get_AIC) {
      Rkk<<-Rk<-list()
      XXtX<- eigenMapMatMult(solve(crossprod(X)),t(X))
      rownames(XXtX)<-colnames(X)
      S=  eigenMapMatMult(X, XXtX)
      for(k in namesX){
        Rkk[[k]]<-Rk[[k]]<-outer(X[,k],XXtX[k,],'*')
      }
      model_lm0@TS<-diag(S)
      AICk<-rep(NA,length(varying))
      names(AICk)<-varying
      TSk<-rep(NA,length(varying))
      TSik<-matrix(0,nrow=n,ncol=length(varying))
      names(TSk)<-colnames(TSik)<-varying
      for(k in varying){
        xtx=eigenMapMatMult(solve(crossprod(X[,k])),t(X[,k]))
        TSik[,k]=diag(eigenMapMatMult(X[,k],xtx))
        AICk[k]<-model_lm0@AIC
      }
      AICc<-AICk[k]
    }
    model_lm0@tS<-sum(model_lm0@TS)
    #browser()

  })
}
