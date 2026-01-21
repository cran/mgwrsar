#' format_and_diagno
#' to be documented
#' @usage format_and_diagno(e)
#' @param e Anenvironment
#' @noRd
#' @return A model of mgwrsar class
format_and_diagno<-function(e){
model=e$model
if(e$TP_estim_as_extrapol) {
  if (!is.null(model$Betav)) XV=e$mymodel@XV[e$TP,]
  if (!is.null(model$Betac)) XC=e$mymodel@XC[e$TP,]
} else {
  if (!is.null(model$Betav)) XV=e$mymodel@XV
  if (!is.null(model$Betac)) XC=e$mymodel@XC
}
names_betav=e$names_betav
names_betav=e$names_betav
new_data=e$new_data
mymodel=e$mymodel
Model=e$Model
isgcv=e$isgcv
family=e$family
pred=e$pred
TP=e$TP
n=e$n
Y=e$Y
isolated_idx=e$isolated_idx
term1 = 0
term2 = 0
if (!is.null(model$Betav) & is.null(new_data) )  term1 <- rowSums(XV * model$Betav)
if (!is.null(model$Betac) & is.null(new_data)) term2 <- XC %*% as.matrix(model$Betac)
if(is.null(new_data) & e$TP_estim_as_extrapol) {residuals <- Y[TP] - term1 - term2
fit = as.numeric(term1 + term2)
} else {residuals <- Y - term1 - term2
fit = as.numeric(term1 + term2)
}
##### DISTINGUER GLM et LM
## CAS GLM
if(Model=='SAR') model$tS=ncol(XC)
if (Model %in% c('GWR_glm','multiscale_gwr') ){
  p<-fit<-family$linkinv(fit)
}
if (!is.null(new_data)){
  if(!(Model %in% c('GWR','MGWR'))){
    if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv'))  {
      lambda_pred=model$Betav[,ncol(model$Betav)]
      if(Model =='MGWRSAR_1_0_kv') beta_pred=model$Betav[,-ncol(model$Betav)] else beta_pred=cbind(matrix(model$Betac,nrow=nrow(model$Betav),ncol=length(model$Betac),byrow=TRUE),model$Betav[,-ncol(model$Betav)])
    } else if(Model %in% c('MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv')) {
      lambda_pred=model$Betac[length(model$Betac)]
      if(Model =='MGWRSAR_0_0_kv')  beta_pred=model$Betav else if(Model=='MGWRSAR_0_kc_kv') beta_pred=cbind(matrix(model$Betac[-length(model$Betac)],nrow=nrow(model$Betav),ncol=length(model$Betac)-1,byrow=TRUE),model$Betav)
    } else if(Model =='MGWRSAR_1_kc_0') {
      lambda_pred=model$Betav[,ncol(model$Betav)]
      beta_pred=matrix(model$Betac,nrow=nrow(model$Betav),ncol=length(model$Betac),byrow=TRUE)
    } else if(Model =='SAR') {
      lambda_pred=model$Betac[length(model$Betac)]
      beta_pred=model$Betac[-length(model$Betac)]
    }
    pred=list(lambda_pred=lambda_pred,beta_pred=beta_pred)
  } else  {
    term1 <- rowSums(XV * model$Betav)
    if(Model=='MGWR') term2 <- new_XC %*% as.matrix(model$Betac)
    pred = as.numeric(term1 + term2)
  }
}
try(colnames(model$Betav) <- names_betav, silent = TRUE)
try(colnames(model$SEV) <- names_betav, silent = TRUE)
try(names(model$Betac) <- names_betac, silent = TRUE)
try(names(model$se) <- names_betac, silent = TRUE)
if(!is.null(new_data)) {mymodel <- list(Betav = model$Betav, Betac = model$Betac,pred = pred,XC = model$XC, XV = model$XV)} else {
  if(!is.null(model$Betav)) mymodel@Betav = model$Betav
  if(!is.null(model$Betac)) mymodel@Betac = model$Betac
  if(!is.null(model$SEV))   mymodel@sev = model$SEV
  if(!is.null(model$se))    mymodel@se = model$se


  ## calcul des residus
  mymodel@isgcv = isgcv
  if(is.null(model$tS)) {
    mymodel@AIC=as.numeric(NA)
    mymodel@AICc=as.numeric(NA)
  } else {
    mymodel@tS    <- model$tS
    if(!is.null(model$Shat)) mymodel@Shat   <- model$Shat
    if(!is.null(model$Rk))  mymodel@R_k  <- model$Rk

    if(Model %in% c('GWR_glm','multiscale_gwr') & family$family %in% c('binomial')){
      SSR=sum(residuals[TP]^2)
      if(Model=='multiscale_gwr') {
        n=n_emp
        mymodel@AIC   <- n*log(SSR/n)+n-model$tS
        mymodel@AICc  <- mymodel@AICc  <- aicc_f(residuals,model$tS,n)
          #n*log(SSR/n)+n*log(2*pi)+n*as.double(n+model$tS)/as.double(n-2-model$tS)+ 2 * model$tS
      } else {
        loglik=sum(Y*log(p) + (1-Y)*log(1-p))
        mymodel@loglik=loglik
        mymodel@AIC   <- -2*loglik+2*(n-model$tS)
        mymodel@AICc  <- mymodel@AICc  <- -2*loglik+2*(model$tS^2+model$tS)/(n-model$tS-2)+ 2 * model$tS
      }
    } else {
      if(is.null(isolated_idx)){
        mymodel@residuals=as.numeric(residuals)
        m=length(TP)
        SSTtp=sum((Y[TP]-mean(Y[TP]))^2)
        SST=sum((Y-mean(Y))^2)
        mymodel@SSRtp<-SSRtp<-sum(mymodel@residuals[TP]^2)
        mymodel@SSR<-SSR<-sum(mymodel@residuals^2)
        mymodel@RMSEtp=sqrt(mean((mymodel@residuals[TP])^2))
        mymodel@RMSE=sqrt(mean((mymodel@residuals)^2)) ## used for optimization with CV
        if(Model!='SAR'){
          if(isgcv) mymodel@CV=mymodel@RMSE
          mymodel@edf   <- n-mymodel@tS*n/m
          mymodel@AIC   <- n*log(SSR/n)+2*(model$tS*n/m)+n+n*log(2*pi)
          mymodel@AICc <- aicc_f(mymodel@residuals,model$tS*n/m,n)
          mymodel@AICctp<- aicc_f(mymodel@residuals[TP],model$tS,m)
          mymodel@R2    <- 1-SSR/SST
          mymodel@R2_adj<- 1-(1-mymodel@R2)*(n-1)/(mymodel@edf-1)
          mymodel@BIC   <- n*log(SSR/n)+n*log(2*pi)+log(n)*model$tS*n/m
          mymodel@loglik <- -(n/2)*(log(SSR/n))-(n/2)*(log(2*pi))-(n/2)
        }
      } else {
        ## if island then use OLS estimate for these obs.
        mymodel@Betav[isolated_idx,]<-matrix(coef(lm.fit(XV,Y)),byrow = T, ncol = ncol(XV),nrow = length(isolated_idx))
        fit[isolated_idx]<-rowSums(mymodel@Betav[isolated_idx,]*XV[isolated_idx,])
        residuals[isolated_idx]<-Y[isolated_idx]-fit[isolated_idx]
        n<-length(TP)
        m<-length(TP)-length(isolated_idx)
        model$tS<-sum(model$TS[-isolated_idx])+length(isolated_idx)/n*ncol(XV)
        mymodel@edf   <- n-mymodel@tS
        mymodel@residuals=as.numeric(residuals)
        mymodel@SSRtp<-SSRtp<-sum(mymodel@residuals[TP]^2)
        mymodel@SSR<-SSR<-sum(mymodel@residuals^2)
        SSTtp=sum((Y[TP]-mean(Y[TP]))^2)
        SST=sum((Y-mean(Y))^2)
        mymodel@RMSEtp<-sqrt(mean((as.numeric(mymodel@residuals)[TP])^2))
        mymodel@RMSE<-sqrt(mean((as.numeric(mymodel@residuals))^2))
        if(isgcv) mymodel@CV=mymodel@RMSE
        m<-n
        mymodel@AIC   <- n*log(SSR/n)+2*(model$tS)+n+n*log(2*pi)
        mymodel@AICctp   <- aicc_f(mymodel@residuals[TP],sum(model$TS[TP][-isolated_idx]),n)## used for optimization with TP
        mymodel@AICc <- aicc_f(mymodel@residuals,model$tS,n)
        mymodel@R2    <- 1-SSR/SST
        mymodel@R2_adj<- 1-(1-mymodel@R2)*(n-1)/(mymodel@edf-1)
        mymodel@BIC   <- n*log(SSR/n)+n*log(2*pi)+log(n)*model$tS*n/m
        mymodel@loglik <- -(n/2)*(log(SSR/n))-(n/2)*(log(2*pi))-(n/2)
      }
    }
  }
  mymodel@fit=fit
  if(!is.null(isolated_idx)) mymodel@isolated_idx=isolated_idx
  my_crs=e$my_crs
  #if(!is.null(pred)) mymodel@pred = pred
}
mymodel
  }
