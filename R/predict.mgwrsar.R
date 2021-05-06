#' mgwrsar Model Predictions
#'
#' predict_mgwrsar is a function for computing predictions of a mgwrsar models. It uses Best Linear Unbiased Predictor for mgwrsar models with spatial autocorrelation.
#' @usage predict_mgwrsar(model, newdata, newdata_coord, W = NULL,
#' type = "BPN", h_w = 100,kernel_w = "knn", k_extra = 12,
#' kernel_extra = "sheppard")
#' @param model   a model of mgwrsar class.
#' @param newdata   a matrix or data.frame of new data.
#' @param newdata_coord   a matrix of new coordinates.
#' @param W   the spatial weight matrix for models with  spatial autocorrelation.
#' @param type  Type for BLUP estimator, default "BPN". If NULL use predictions without spatial bias correction.
#' @param h_w   bandwidth for constructing W, if W is NULL.
#' @param kernel_w   kernel for constructing W, if W is NULL.
#' @param k_extra   number of neighboors for local parameter extrapolation, default 12.
#' @param kernel_extra kernel for local parameter extrapolation, default sheppard kernel.
#'
#' @return A vector of predictions.
#' @seealso  MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, kernelW_C
#' @examples
#' \donttest{
#' library(mgwrsar)
#' data(mydata)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' W=KNN(coord,2)
#' model_GWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[1:800,],
#' coord=coord[1:800,],fixed_vars=NULL,kernels=c('gauss_adapt'),H=50,
#'  Model = 'GWR',control=list())
#' Y_pred=predict_mgwrsar(model_GWR_insample, newdata=mydata[801:1000,],
#' newdata_coord=coord[801:1000,],k_extra = 8, kernel_extra = "sheppard")
#' head(Y_pred)
#' head(mydata$Y_gwr[801:1000])
#' sqrt(mean((mydata$Y_gwr[801:1000]-Y_pred)^2))
#'
#' ## predict with spatial autocorrelation
#' model_MGWRSAR_1_0_kv_insample<-MGWRSAR(formula = 'Y_mgwrsar_1_0_kv~X1+X2+X3',data = mydata[1:800,],
#' coord=coord[1:800 ,], fixed_vars=NULL,kernels=c('gauss_adapt'),H=50,
#' Model = 'MGWRSAR_1_0_kv',control=list(W=W[1:800,1:800],Lambdacor=TRUE,SE=TRUE))
#' summary_mgwrsar(model_MGWRSAR_1_0_kv_insample)
#'
#' ## with BLUP
#' Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=mydata[801:1000,],
#' newdata_coord=coord[801:1000,], k_extra = 12,  W = W,
#'  type = "BPN",kernel_extra = "sheppard")
#' head(Y_pred)
#' head(mydata$Y_gwr[801:1000])
#' sqrt(mean((mydata$Y_gwr[801:1000]-Y_pred)^2))
#'
#' ## without BLUP
#' Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=mydata[801:1000,],
#' newdata_coord=coord[801:1000,], k_extra = 12,  W = W,
#' type = "TC",kernel_extra = "sheppard")
#' head(Y_pred)
#' head(mydata$Y_mgwrsar_1_0_kv[801:1000])
#' sqrt(mean((mydata$Y_mgwrsar_1_0_kv[801:1000]-Y_pred)^2))
#' }
predict_mgwrsar  <- function(model, newdata, newdata_coord, W = NULL, type = "BPN", h_w = 100,kernel_w = "knn", k_extra = 12, kernel_extra = "sheppard") {
if(class(model)!='mgwrsar') stop('A mgwrsar class object is needed')
if(is.null(newdata) | is.null(newdata_coord)) stop('provide newdata and newdata_coord objects')
#######
insample_size=length(model$Y) ## insample
outsample_size=nrow(newdata)  ###outsample
if(ncol(newdata_coord)>2) coords=rbind(model$S,newdata_coord) else coords=rbind(model$coord,newdata_coord)
m=insample_size
n=outsample_size+insample_size
S=1:m
O=(m+1):(n)

if(!('Intercept' %in% colnames(newdata))) {
newdata=cbind(rep(1,nrow(newdata)),newdata)
colnames(newdata)[1]='Intercept'
}

###
if(model$Model=='MGWR') {
  newdata=newdata[,c(colnames(model$Betav),names(model$Betac))]
  } else if(model$Model=='MGWRSAR_0_0_kv') {
  newdata=newdata[,c(colnames(model$Betav),names(model$Betac)[-length(model$Betac)])]
  beta_hat=model$Betav
  lambda_hat=model$Betac
  } else if(model$Model == 'MGWRSAR_1_0_kv') {
    newdata=newdata[,c(colnames(model$Betav)[-ncol(model$Betav)],names(model$Betac))]
    beta_hat=model$Betav[,-ncol(model$Betav)]
    lambda_hat=model$Betav[,ncol(model$Betav)]
  }
if(model$Model=='SAR') {
  newdata=newdata[,names(model$Betac)[-length(names(model$Betac))]]
  beta_hat=model$Betac[-length(model$Betac)]
  lambda_hat=model$Betac[length(model$Betac)]
}

if(!is.null(model$h_w)){
h_w=model$h_w
kernel_w=model$kernel_w
}

if(is.null(W)) {
  if(ncol(coords)==2) W<-KNNX(coords,k_neighbors=h_w,diagnull=TRUE,kernel=kernel_w,query=NULL) else {
    stop( "To be coded")
  }
}
YS = model$Y
e=model$residuals
###  X = rbind(X_insample,X_outsample)
###  e=residuals_insample
###  beta_hat,lambda_hat = coef_insample
###  S,O index of insample and outsample
###  dim(W)=length(c(S,O)) x length(c(S,O))


if(model$Model=='SAR') Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model$X,newdata)),W,e,beta_hat=beta_hat,lambda_hat=lambda_hat,S,O,type,coords=coords)

if(model$Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv'))  Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model$X,newdata)),W,e,beta_hat=beta_hat,lambda_hat=lambda_hat,S,O,type,k_extra=k_extra,kernel_extra=kernel_extra,model='MGWRSARC',coords=coords)

if(model$Model=='GWR') {
  X=as.matrix(rbind(model$X,setNames(cbind(rep(1,nrow(newdata)),newdata[,colnames(newdata)  %in% colnames(model$X)[-1]]),colnames(model$X))))
W_extra=KNNX(coord[S,],k_neighbors=k_extra,query=coords[O,],kernel=kernel_extra)

Beta_proj_out=as.matrix(W_extra %*% model$Betav)
Betav=rbind(model$Betav,Beta_proj_out)
Y_predicted=rowSums(as.matrix(Betav*X))
Y_predicted=Y_predicted[O]
}

if(model$Model=='MGWR') {
  X=as.matrix(rbind(model$X,setNames(cbind(rep(1,nrow(newdata)),newdata[,colnames(newdata)  %in% colnames(model$X)[-1]]),colnames(model$X))))
W_extra=KNNX(coord[S,],k_neighbors=k_extra,query=coords[O,],kernel=kernel_extra)
betav=as.matrix(cbind(model$Betav,matrix(model$Betac,nrow=m,ncol=length(model$Betac), byrow =T)))
Beta_proj_out=as.matrix(Wk %*% betav)
Betav=rbind(betav,Beta_proj_out)
Y_predicted=rowSums(as.matrix(Betav*X))
Y_predicted=Y_predicted[O]
}
Y_predicted
}
