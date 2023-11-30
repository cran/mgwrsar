#' mgwrsar Model Predictions
#' predict_mgwrsar is a function for computing predictions of a mgwrsar models. It uses Best Linear Unbiased Predictor for mgwrsar models with spatial autocorrelation.
#' @usage predict_mgwrsar(model, newdata, newdata_coords, W = NULL, type = "BPN",
#' h_w = 100,kernel_w = "rectangle",maxobs=4000,beta_proj=FALSE,
#' method_pred='TP', k_extra = 8)
#' @param model a model of mgwrsar class.
#' @param newdata a matrix or data.frame of new data.
#' @param newdata_coords  a matrix of new coordinates, and eventually other variables if a General Kernel Product is used.
#' @param W the spatial weight matrix for models with  spatial autocorrelation.
#' @param type Type for BLUP estimator, default "BPN". If NULL use predictions without spatial bias correction.
#' @param  h_w A bandwidth value for the spatial weight matrix
#' @param kernel_w kernel type for the spatial weight matrix. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"),
#' epanechnikov ("epane"), gaussian ("gauss")) .
#' @param maxobs  maximum number of observations for exact calculation of solve(I- rho*W), default maxobs=4000.
#' @param beta_proj A boolean, if TRUE the function then return a two elements list(Y_predicted,Beta_proj_out)
#' @param method_pred If method_pred = 'TP' (default) prediction is done by recomputing a MGWRSAR model
#' with new-data as target points, else if method_pred in ('tWtp_model','model','sheppard') a matrix
#' for projecting estimated betas is used (see details).
#' @param k_extra number of neighboors for local parameter extrapolation if sheppard kernel is used, default 8.
#' @return A vector of predictions if beta_proj is FALSE or a list with a vector named Y_predicted and a matrix named Beta_proj_out.
#' @details if method_pred ='tWtp_model',  the weighting matrix for prediction is
#' based on the expected weights of outsample data if they were had been added to
#' insample data to estimate the corresponding MGWRSAR (see Geniaux 2022 for
#' further detail), if method_pred ='sheppard'a sheppard kernel with k_extra neighbours (default 8) is used and if method_pred='kernel_model' the same kernel
#' and number of neighbors as for computing the MGWRSAR model is used.
#' @seealso  MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, kernel_matW
#' @examples
#' \donttest{
#' library(mgwrsar)
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#' length_out=800
#' index_in=sample(1:1000,length_out)
#' index_out=(1:1000)[-index_in]
#'
#' model_GWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[index_in,],
#' coords=coords[index_in,],fixed_vars=NULL,kernels=c ('gauss'),H=8, Model = 'GWR',
#' control=list(adaptive=TRUE))
#' summary_mgwrsar(model_GWR_insample)
#'
#' newdata=mydata[index_out,]
#' newdata_coords=coords[index_out,]
#' newdata$Y_mgwrsar_1_0_kv=0
#'
#' Y_pred=predict_mgwrsar(model_GWR_insample, newdata=newdata,
#' newdata_coords=newdata_coords)
#' head(Y_pred)
#' head(mydata$Y_gwr[index_out])
#' sqrt(mean((mydata$Y_gwr[index_out]-Y_pred)^2)) # RMSE
#' }
predict_mgwrsar  <- function(model, newdata, newdata_coords, W = NULL, type = "BPN", h_w = 100,kernel_w = "rectangle",maxobs=4000,beta_proj=FALSE,method_pred='TP', k_extra = 8) {
  if(method_pred !='sheppard' & model$Model %in% c('multiscale_gwr','tds_mgwr')) {
    cat("Warning: only sheppard method for method_pred can be used with multiscale_gwr and tds_mgwr \n automatic switch to method_pred='sheppard'")
    method_pred='sheppard'
  }
  if(is.null(model$TP) | length(model$TP)==nrow(model$X)) noTP=TRUE else noTP=FALSE
if(!is(model,'mgwrsar')) stop('A mgwrsar class object is required')
if(model$Model %in% c('MGWR','MGWRSAR_1_0_kv','MGWRSAR_1_kc_0','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv') & method_pred=='TP') {
  cat("Warnings: method_pred=='TP' is not implemented for  Model in ('MGWR','MGWRSAR_1_kc_0','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv'): automatic swith to method_pred='tWtp_model'")
  method_pred='tWtp_model'
}
if(length(model$TP)!=length(model$Y) & method_pred=='TP') {
    cat("Warnings: if estimated model used target points, method_pred ='TP' may be innapropriate if out-sample size is large: method_pred='tWtp_model' should be faster")
}

if(is.null(newdata) | is.null(newdata_coords)) stop('provide the newdata and newdata_coords objects')
#######
  newdata<-data.frame(newdata)
  myYname=terms(as.formula(model$formula))[[2]]
  ## add a null variable with the response variable name
  if(!(as.character(myYname) %in% colnames(newdata))){
  newdata[[myYname]]<-0
}
####
  insample_size=length(model$Y) ## insample
  outsample_size=nrow(newdata)  ###outsample

  if(ncol(newdata_coords)>2) {
    colnames(model$S)<- colnames(newdata_coords)
    coords=rbind(model$S,newdata_coords)
    } else {
      colnames(model$coords)<- colnames(newdata_coords)
      coords=rbind(model$coords,newdata_coords)
    }
  m=insample_size
  n=outsample_size+insample_size
  S=1:m
  O=(m+1):(n)
  YS = model$Y
  e=model$residuals

  ## Use all coefficients instead of coefficients located in TP seems preferable. To be checked with real data
  #model$TP=1:length(model$Y)

  if(method_pred=='TP' & (model$Model %in% c('GWR'))){    #,'MGWRSAR_1_0_kv'
  TP=NULL
  #if(model$Model =='GWR') {
    model$mycall$control[['W']]<-quote(model$W)
    con=eval(model$mycall$control)
    #con=c(con, list(S_out=newdata_coords,new_data=newdata,new_W=W))
    full_data=rbind(newdata,model$data)
    full_coords=rbind(newdata_coords,model$coords)
    TP=1:length(O)
    con=c(con, list(TP=TP,isgcv=TRUE,S_out=TRUE))
    pred_TP<-MGWRSAR(formula = model$formula, data = full_data,coords=full_coords, fixed_vars=model$fixed_vars,kernels=model$kernels,H=model$H, Model=model$Model ,control=con)
    Y_predicted=pred_TP$fit[TP]
    # } else if(model$Model =='MGWRSAR_1_0_kv'){
    # lambda_in=model$Betav[,ncol(model$Betav)]
    # beta_in=model$Betav[,-ncol(model$Betav)]
    # model$mycall$control[['W']]<-quote(model$W)
    # con=call_modify(model$mycall$control, S_out=quote(newdata_coords),new_data=quote(newdata),new_W=quote(W))
    # pred_TP<-MGWRSAR(formula = model$formula, data = model$data,coords=model$coords, fixed_vars=model$fixed_vars,kernels=model$kernels,H=model$H, Model=model$Model ,control=eval(con))
    #
    # beta_out=pred_TP$pred$beta_pred
    # lambda_out=pred_TP$pred$lambda_pred
    # X=rbind(model$XV[,-ncol(model$XV)],pred_TP$XV)
    # Y_predicted=BP_pred_SAR(YS=model$Y, X=X, W=W, e=model$residuals, beta_hat=rbind(beta_in,beta_out), lambda_hat=c(lambda_in,lambda_out), S, O, type = type,model =model$Model)
    # Betav_proj_out=pred_TP$Betav
    # Betac_proj_out=pred_TP$Betac
    # }
  } else
    {
#######
if(!('Intercept' %in% colnames(newdata))) {
newdata=cbind(rep(1,nrow(newdata)),newdata)
colnames(newdata)[1]='Intercept'
}

if(!(model$Model %in% c('OLS','SAR'))){
  ###
  if(model$Model %in% c('GWR','multiscale_gwr','tds_mgwr')) {
    beta_in=model$Betav
  } else if(model$Model=='MGWR') {
    beta_in=cbind(model$Betav,matrix(model$Betac,nrow=m,ncol=length(model$Betac),byrow=TRUE))
    colnames(beta_in)=c(colnames(model$Betav),names(model$Betac))
  } else if(model$Model=='MGWRSAR_1_0_kv') {
    lambda_in=model$Betav[,ncol(model$Betav)]
    beta_in=model$Betav[,-ncol(model$Betav)]
  } else if(model$Model=='MGWRSAR_1_kc_kv') {
    lambda_in=model$Betav[,ncol(model$Betav)]
    beta_in=cbind(model$Betav[,-ncol(model$Betav)],matrix(model$Betac,nrow=m,ncol=length(model$Betac),byrow=TRUE))
    colnames(beta_in)=c(colnames(model$Betav)[-ncol(model$Betav)],names(model$Betac))
  } else if(model$Model=='MGWRSAR_1_kc_0') {
    lambda_in=model$Betav
    beta_in=matrix(as.matrix(model$Betac),nrow=m,ncol=length(model$Betac),byrow=TRUE)
    colnames(beta_in)=names(model$Betac)
  } else if(model$Model %in% c('MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv')){
    lambda_in=rep(model$Betac[length(model$Betac)],m)
    if(model$Model=='MGWRSAR_0_0_kv') beta_in=model$Betav  else if(model$Model=='MGWRSAR_0_kc_kv') {
      beta_in=cbind(model$Betav,matrix(model$Betac[-length(model$Betac)],nrow=m,ncol=length(model$Betac)-1,byrow=TRUE))
      colnames(beta_in)=c(colnames(model$Betav),names(model$Betac)[-length(model$Betac)])
    }
  }
  #beta_in<-beta_in[model$TP,colnames(model$X)]
  namesX=colnames(model$X)
  if(!(model$Model %in% c('GWR','MGWR','multiscale_gwr','tds_mgwr'))){
    lambda_in<-lambda_in[model$TP]
    il=which(namesX=='lambda')
    if(length(il)>0) namesX<-namesX[-il]
    }
  newdata=newdata[,namesX]
  #if(model$Model =='multiscale_GWR') newdata=apply(newdata,2,function(x) scale(x,scale=FALSE))

  if(method_pred=='tWtp_model'){
    W_extra=kernel_matW(H=model$H,kernels=model$kernels,coord_i=coords[S,],coord_j=coords,NN=model$NN,Type=model$Type,adaptive=model$adaptive,diagnull=FALSE,rowNorm=TRUE)[,-(1:nrow(model$Betav))]
    W_extra<- normW(Matrix::t(W_extra))
  } else if(method_pred=='model'){
    W_extra=kernel_matW(H=model$H,kernels=model$kernels,coord_i=rbind(model$coords,newdata_coords),coord_j=model$coords,NN=model$NN,Type=model$Type,adaptive=model$adaptive,diagnull=FALSE,rowNorm=TRUE)[-(1:nrow(model$Betav)),]
  } else {
    W_extra=kernel_matW(H=k_extra,kernels='sheppard',coord_i=rbind(model$coords,newdata_coords),coord_j=model$coords,NN=k_extra+1,Type=model$Type,adaptive=FALSE,diagnull=FALSE,rowNorm=TRUE)[-(1:nrow(model$Betav)),]

  }
} else if(model$Model=='SAR') {
  newdata=newdata[,names(model$Betac)[-length(names(model$Betac))]]
  beta_in=model$Betac[-length(model$Betac)]
  lambda_in=model$Betac[length(model$Betac)]
} else if(model$Model=='OLS'){
  newdata=newdata[,names(model$Betac)]
  Y_predicted=as.matrix(newdata) %*% model$Betac
}
if(!is.null(model$h_w)){
h_w=model$h_w
kernel_w=model$kernel_w
}

if(is.null(W) & !(model$Model %in% c('OLS','GWR','MGWR','multiscale_gwr','tds_mgwr'))) {
  W=kernel_matW(H=h_w,kernels=kernel_w,coord_i=coords,diagnull=TRUE,NN=h_w,adaptive=TRUE,rowNorm=TRUE)
}

if(model$Model=='SAR') {
  Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model$X,newdata)),W,e,beta_hat=beta_in,lambda_hat=lambda_in,S,O,type,coords=coords,maxobs=maxobs)

  } else if(model$Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_1_kc_0'))  {
    Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model$X,newdata)),W,e,beta_hat=rbind(beta_in,as.matrix(W_extra %*%beta_in)),lambda_hat=c(lambda_in,as.numeric(W_extra %*%lambda_in)),S,O,type,k_extra=k_extra,kernel_extra=kernel_extra,model=model$Model, W_extra = W_extra,coords=coords,maxobs=maxobs)

  } else if (model$Model %in% c('GWR','GWRtp','GWRboost','MGWR','multiscale_gwr','tds_mgwr')) {
    X=as.matrix(newdata)
    Beta_proj_out=as.matrix(W_extra %*% beta_in)
    #else Beta_proj_out=as.matrix(W_extra %*% beta_in[model$TP,])
    #Betav=rbind(beta_in,Beta_proj_out)
    Y_predicted=rowSums(as.matrix(Beta_proj_out*X))
    #Y_predicted=Y_predicted[O]
  }
}
if(beta_proj) return(list(Y_predicted=Y_predicted,Betav_proj_out=Betav_proj_out,Betac_proj_out=Betac_proj_out))
else return(Y_predicted)
}
