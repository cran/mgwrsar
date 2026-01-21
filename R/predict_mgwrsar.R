#' mgwrsar Model Predictions
#' predict_mgwrsar is a function for computing predictions of a mgwrsar models. It uses Best Linear Unbiased Predictor for mgwrsar models with spatial autocorrelation.
#' @usage predict_mgwrsar(model, newdata, newdata_coords, W = NULL,
#' type = "BPN", h_w = 100,kernel_w = "rectangle",maxobs=4000,
#' beta_proj=FALSE,method_pred='TP', k_extra = 8,exposant=6)
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
#' with new-data as target points, else if method_pred in ('tWtp_model','model','shepard') a matrix
#' for projecting estimated betas is used (see details).
#' @param k_extra number of neighboors for local parameter extrapolation if shepard kernel is used, default 8.
#' @param exposant shapenig parameter for tds_mgtwr model, default 6.
#' @return A vector of predictions if beta_proj is FALSE or a list with a vector named Y_predicted and a matrix named Beta_proj_out.
#' @details if method_pred ='tWtp_model',  the weighting matrix for prediction is
#' based on the expected weights of outsample data if they were had been added to
#' insample data to estimate the corresponding MGWRSAR (see Geniaux 2022 for
#' further detail), if method_pred ='shepard'a shepard kernel with k_extra neighbours (default 8) is used and if method_pred='kernel_model' the same kernel
#' and number of neighbors as for computing the MGWRSAR model is used.
#' @seealso  MGWRSAR, golden_search_bandwidth, plot_mgwrsar, kernel_matW
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
#' summary(model_GWR_insample)
#'
#' newdata=mydata[index_out,]
#' newdata_coords=coords[index_out,]
#' newdata$Y_mgwrsar_1_0_kv=0
#'
#' Y_pred=predict(model_GWR_insample, newdata=newdata,
#' newdata_coords=newdata_coords)
#' head(Y_pred)
#' head(mydata$Y_gwr[index_out])
#' sqrt(mean((mydata$Y_gwr[index_out]-Y_pred)^2)) # RMSE
#' }
#' @noRd
predict_mgwrsar  <- function(model, newdata, newdata_coords, W = NULL, type = "BPN", h_w = 100,kernel_w = "rectangle",maxobs=4000,beta_proj=FALSE,method_pred='TP', k_extra = 8,exposant=8) {
  ##############
  ##### WARNINGS
  ##############
  # if(!(method_pred %in% c('Multiscale_GDT','shepard','tWtp_model')) & model@Model %in% c('multiscale_gwr','tds_mgwr','atds_gwr','atds_mgwr')) {
  #   cat("Warning: only shepard and tWtp_model method for method_pred can be used with multiscale_gwr and tds_mgwr \n automatic switch to method_pred='tWtp_model'")
  #   method_pred='tWtp_model'
  # }

  if(is.null(model@TP) | length(model@TP)==nrow(model@X)) noTP=TRUE else noTP=FALSE

if(!noTP & method_pred=='TP') {
  cat("Warnings: method_pred=='TP' is not recommanded when the model has already been estimated using TP,  automatic swith to method_pred='tWtp_model'")
  method_pred='tWtp_model'
}


if(!is(model,'mgwrsar')) stop('A mgwrsar class object is required')

  tp_supported <- model@Model %in% c("GWR","MGWR")  # + autres si tu les supportes vraiment

  if (method_pred == "TP" && !tp_supported) {
    warning("method_pred='TP' is only implemented for ... Switching to 'tWtp_model'.")
    method_pred <- "tWtp_model"
  }

# if(model@Model %in% c('MGWR','MGWRSAR_1_0_kv','MGWRSAR_1_kc_0','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv') & method_pred=='TP') {
#   cat("Warnings: method_pred=='TP' is not implemented for  Model in ('MGWR','MGWRSAR_1_kc_0','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv'): automatic swith to method_pred='tWtp_model'")
#   method_pred='tWtp_model'
# }

if(model@Model %in% c("multiscale_gwr", "tds_mgwr", "atds_mgwr", "tds_gwr") & method_pred=='TP') {
    cat("Warnings: method_pred=='TP' is not implemented for  Model in (multiscale_gwr, tds_mgwr, atds_mgwr, tds_gwr): automatic swith to method_pred='shepard'")
    method_pred='shepard'
  }
if(length(model@TP)!=length(model@Y) & method_pred=='TP') {
    cat("Warnings: if estimated model used target points, method_pred ='TP' may be innapropriate if out-sample size is large: method_pred='tWtp_model' should be faster")
}
if(model@Type=='GDT' &  !(method_pred %in% c('model','tWtp_model'))) stop('Only method="model" is available for prediction with space-time kernel (Type="GDT")')

if(is.null(newdata) | is.null(newdata_coords)) stop('provide the newdata and newdata_coords objects')

  ##############
  ##### O S prep
  ##############

  if(!(any(c('Intercept','(Intercept)') %in% colnames(newdata)))) {
    newdata=cbind(rep(1,nrow(newdata)),newdata)
    colnames(newdata)[1]='Intercept'
  } else if(colnames(newdata)[1]=='(Intercept)') colnames(newdata)[1]='Intercept'
  names_newX=colnames(newdata)
  newdata<-data.frame(newdata)
  colnames(newdata)<-names_newX
  myYname=terms(as.formula(model@formula))[[2]]
  if(!(as.character(myYname) %in% colnames(newdata))){
  newdata[[myYname]]<-0
}
  insample_size=length(model@Y)
  outsample_size=nrow(newdata)


  if(model@Type=='GDT') {

      if(sum(duplicated(newdata_coords[,1:2]))>0) {
        newdata_coords[,1] <- newdata_coords[,1] + seq_along(newdata_coords[,1]) * 1e-12
        newdata_coords[,2] <- newdata_coords[,2] + seq_along(newdata_coords[,2]) * 1e-12
      }
    if(sum(duplicated(newdata_coords[,3]))>0) {
      newdata_coords[,3] <- newdata_coords[,3] + seq_along(newdata_coords[,3]) * 1e-12
    }
    A <- as.data.frame(cbind(model@coords, model@Z))
    names(A) <- c("x","y","time")
    B <- as.data.frame(newdata_coords)
    names(B) <- c("x","y","time")
    A$time <- as.numeric(A$time)
    B$time <- as.numeric(B$time)
    coords <- rbind(A, B)
    } else if(model@Type=='T') {
      coords=as.matrix(c(model@Z,as.matrix(newdata_coords)[,ncol(newdata_coords)]),ncol=1)
    } else {
      colnames(model@coords)<- colnames(newdata_coords)[1:2]
      coords=rbind(model@coords,newdata_coords[,1:2])
    }
  coords<-make_unique_by_structure(coords)
  m=insample_size
  n=outsample_size+insample_size
  S=1:m
  O=(m+1):(n)
  YS = model@Y
  e=model@residuals
  namesX=colnames(model@X)


  if(is.null(W) & !(model@Model %in% c('OLS','GWR','GWR_glmboost','MGWR','multiscale_gwr','tds_mgwr','atds_mgwr','atds_gwr'))) {
    if(!is.null(model@h_w)){
      h_w=model@h_w
      kernel_w=model@kernel_w
    }
    W=kernel_matW(H=h_w,kernels=kernel_w,coords=coords,diagnull=TRUE,NN=h_w,adaptive=TRUE)
  }

  ##############
  ##### ESTIMATION OR EXTAPOLATION
  ##############

  if(method_pred=='TP' & (model@Model %in% c('GWR','MGWR'))){    #,'MGWRSAR_1_0_kv'

    ##############
    ##### ESTIMATION
    ##############
  TP=NULL
  #if(model@Model =='GWR') {
  # ============================================
  # FIX : Always ensure mycall$control is a list
  # ============================================
  ctrl <- model@mycall$control
  if (is.null(ctrl) || !is.list(ctrl)) {
    ctrl <- list()
  }
  ctrl$W <- quote(model@W)
  model@mycall$control <- ctrl
    model@mycall$control[['W']]<-quote(model@W)
    con=eval(model@mycall$control)
    model@data$Intercept=1
    full_data=rbind(newdata[,names(model@data)],model@data)
    TP=1:length(O)
    con=c(con, list(TP=TP,isgcv=TRUE,TP_estim_as_extrapol=TRUE,get_ts=FALSE,NN=model@NN))
    if(ncol(newdata_coords)>2) {
      full_coords=rbind(newdata_coords[,-3],model@coords)
      con[['Z']]=c(newdata_coords[,3],model@Z)
      } else full_coords=rbind(newdata_coords,model@coords)
    pred_TP<-MGWRSAR(formula = model@formula, data = full_data,coords=full_coords[,-3], fixed_vars=model@fixed_vars,kernels=model@kernels,H=model@H, Model=model@Model ,control=con)
    Beta_proj_out=pred_TP@Betav
    Y_predicted=pred_TP@fit[TP]
    # } else if(model@Model =='MGWRSAR_1_0_kv'){
    # lambda_in=model@Betav[,ncol(model@Betav)]
    # beta_in=model@Betav[,-ncol(model@Betav)]
    # model@mycall$control[['W']]<-quote(model@W)
    # con=call_modify(model@mycall$control, TP_estim_as_extrapol=quote(newdata_coords),new_data=quote(newdata),new_W=quote(W))
    # pred_TP<-MGWRSAR(formula = model@formula, data = model@data,coords=model@coords, fixed_vars=model@fixed_vars,kernels=model@kernels,H=model@H, Model=model@Model ,control=eval(con))
    #
    # beta_out=pred_TP$pred$beta_pred
    # lambda_out=pred_TP$pred$lambda_pred
    # X=rbind(model@XV[,-ncol(model@XV)],pred_TP$XV)
    # Y_predicted=BP_pred_SAR(YS=model@Y, X=X, W=W, e=model@residuals, beta_hat=rbind(beta_in,beta_out), lambda_hat=c(lambda_in,lambda_out), S, O, type = type,model =model@Model)
    # Betav_proj_out=pred_TP$Betav
    # Betac_proj_out=pred_TP$Betac
    # }
  } else {
##############
##### EXTAPOLATION
##############

##############
####### beta_in / lambda_in prep
##############
if(model@Model=='OLS') {
  newdata=newdata[,names(model@Betac)]
  Y_predicted=as.matrix(newdata) %*% model@Betac
  Beta_proj_out=matrix(model@Betac,ncol=length(model@Betac),nrow=nrow(newdata),byrow=TRUE)

} else if(model@Model=='SAR'){
  newdata=newdata[,names(model@Betac)[-length(names(model@Betac))]]
  newdata=newdata[,names(model@Betac)[-length(names(model@Betac))]]
  beta_in=model@Betac[-length(model@Betac)]
  lambda_in=model@Betac[length(model@Betac)]
  Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model@X,newdata)),W,e,beta_hat=beta_in,lambda_hat=lambda_in,S,O,type,coords=coords,maxobs=maxobs)
  Beta_proj_out=matrix(beta_in,ncol=length(beta_in),nrow=nrow(newdata),byrow=TRUE)


} else if(model@Model %in% c('GWR',"GWR_glm","GWR_glmboost",'GWR_gamboost_linearized')){
  beta_in=model@Betav

} else if(model@Model=='MGWR'){
  beta_in=cbind(model@Betav,matrix(model@Betac,nrow=m,ncol=length(model@Betac),byrow=TRUE))
  colnames(beta_in)=c(colnames(model@Betav),names(model@Betac))

} else if(model@Model %in% c('multiscale_gwr','tds_mgwr','atds_mgwr','atds_gwr')){
beta_in=model@Betav

} else if(model@Model=='MGWRSAR_0_0_kv'){
  beta_in=model@Betav
  lambda_in=rep(model@Betac[length(model@Betac)],m)
} else if(model@Model=='MGWRSAR_0_kc_kv'){
  lambda_in=rep(model@Betac[length(model@Betac)],m)
  beta_in=cbind(model@Betav,matrix(model@Betac[-length(model@Betac)],nrow=m,ncol=length(model@Betac)-1,byrow=TRUE))
  colnames(beta_in)=c(colnames(model@Betav),names(model@Betac)[-length(model@Betac)])
} else if(model@Model=='MGWRSAR_1_kc_kv'){
  lambda_in=model@Betav[,ncol(model@Betav)]
  beta_in=cbind(model@Betav[,-ncol(model@Betav)],matrix(model@Betac,nrow=m,ncol=length(model@Betac),byrow=TRUE))
  colnames(beta_in)=c(colnames(model@Betav)[-ncol(model@Betav)],names(model@Betac))
} else if(model@Model=='MGWRSAR_1_0_kv'){
  lambda_in=model@Betav[,ncol(model@Betav)]
  beta_in=model@Betav[,-ncol(model@Betav)]
} else if(model@Model=='MGWRSAR_1_kc_0'){
  lambda_in=model@Betav
  beta_in=matrix(as.matrix(model@Betac),nrow=m,ncol=length(model@Betac),byrow=TRUE)
  colnames(beta_in)=names(model@Betac)
}
if(!(model@Model %in% c('OLS','SAR'))) tokeep=which(!is.na(beta_in[,1]))

##############
##### calcul W_extra et Beta_proj_out
##############

if(model@Model %in% c('multiscale_gwr','tds_mgwr','atds_mgwr','atds_gwr')){
  if(method_pred!='tWtp_model') res=prep_d(coords,NN=model@NN,TP=O,extrapol=TRUE,ratio=1,QP=S,kernels=model@kernels,Type=model@Type)
  if(method_pred!='shepard'){
  Beta_proj_out=matrix(0,ncol=ncol(beta_in),nrow=length(O))
  colnames(Beta_proj_out)=colnames(beta_in)
  for(k in colnames(beta_in)){ #parallel ?
    #cat(k,' ')
      if(method_pred=='tWtp_model'){
        #browser()
        if(length(model@Ht)>0) myH=c(model@H[k],model@Ht[k]) else myH=model@H[k]
        W_extra=kernel_matW(H=myH,kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=S,Type=model@Type,dists=res$dists,diagnull=FALSE,extrapol=F,alpha=model@alpha)[,O]
        # W_extra = t(apply(W_extra, 1, function(x) x * as.numeric(x >= sort(x, decreasing = T)[k_extra + 2])))
        W_extra<-W_extra[tokeep,]
      } else if(method_pred=='model'){
        if(length(model@Ht)>0) {
          myH=c(model@H[k],model@Ht[k])
          if(all(myH==c(model@V[1],model@Vt[1]))) {
            W_extra=Matrix(1,nrow=length(O),ncol=length(S))
          } else if(myH[2]==model@Vt[1]) {
            if(!exists('resS')){
              resS=prep_d(coords[c(S,O),1:2],NN=model@NN,TP=O,extrapol=TRUE,ratio=1,QP=S,kernels=model@kernels[1],Type='GD')
            }
            W_extra=kernel_matW(H=myH[1],kernels=model@kernels[1],coords=coords[c(S,O),1:2],NN=model@NN,TP=O,Type='GD',adaptive=model@adaptive[1],indexG=resS$indexG,dists=resS$dists,diagnull=FALSE,extrapol=T,QP=S,alpha=model@alpha)
          } else if(myH[1]==model@V[1]) {
            if(!exists('resT')){
              resT=prep_d(as.matrix(coords[c(S,O),3],ncol=1),NN=model@NN,TP=O,extrapol=TRUE,ratio=1,QP=S,kernels=model@kernels[2],Type='T')
            }
            W_extra=kernel_matW(H=myH[2],kernels=model@kernels[2],coords=as.matrix(coords[c(S,O),3],ncol=2),NN=model@NN,TP=O,Type='T',adaptive=F,indexG=resT$indexG,dists=resT$dists,diagnull=FALSE,extrapol=T,QP=S,alpha=model@alpha)
            } else {

              W_extra=kernel_matW(H=myH,kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,indexG=res$indexG,dists=res$dists,diagnull=FALSE,extrapol=T,QP=S,alpha=model@alpha)
            }
         # W_extra=kernel_matW(H=myH,kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,indexG=res$indexG,dists=res$dists,diagnull=FALSE,extrapol=T,QP=S,alpha=model@alpha) # visiblement c'est mieux en pure GDT pour les prédictions
        } else {
          myH=model@H[k]
          if(myH[1]==model@V[1])  W_extra=matrix(1,nrow=length(O),ncol=length(S)) else W_extra=kernel_matW(H=myH,kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,indexG=res$indexG,dists=res$dists,diagnull=FALSE,extrapol=T,QP=S,alpha=model@alpha)
        }
        #if(model@Type=='GDT') W_extra = t(apply(W_extra, 1, function(x) x * as.numeric(x >= sort(x, decreasing = T)[k_extra + 2])))
        W_extra<-normW(W_extra^exposant)
      }
      Beta_proj_out[,k]=as.matrix(W_extra %*% beta_in[tokeep,k])
  }
  if(exists('res')) rm(res)
  if(exists('resS')) rm(resS)
  if(exists('resT')) rm(resT)
  gc()
  } else if(method_pred=='shepard'){
    if(model@Type=='GD'){
      W_extra=kernel_matW(H=k_extra,kernels='shepard',coords=coords[c(S,O),],NN=k_extra+1,TP=O,Type=model@Type,adaptive=F,diagnull=FALSE,extrapol=T,QP=S)
    } else {
      range01 <- function(x){(x-min(x))/(max(x)-min(x))};
      coords2=apply(coords,2,range01)
      W_extra=kernel_matW(H=k_extra,kernels='shepard',coords=coords2[c(S,O),],NN=k_extra+1,TP=O,Type='GD',adaptive=F,diagnull=FALSE,extrapol=T,QP=S)
    }
    Beta_proj_out=as.matrix(W_extra %*% beta_in[tokeep,])
  }
} else if(model@Model %in% c('GWR','MGWR','MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv','MGWRSAR_1_kc_kv','MGWRSAR_1_0_kv','MGWRSAR_1_kc_0',"GWR_glm","GWR_glmboost",'GWR_gamboost_linearized')){

  if(method_pred=='tWtp_model'){
    if(model@Type=='GD') {
      W_extra=kernel_matW(H=c(model@H,model@Ht),kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=S,Type=model@Type,adaptive=model@adaptive,diagnull=FALSE,extrapol=F)[,O]
      W_extra<-W_extra[tokeep,]
      W_extra<- normW(Matrix::t(W_extra))
    } else if(model@Type=='GDT'){
        W_extra=kernel_matW(H=c(model@H,model@Ht),kernels=model@kernels,coords=coords,NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,diagnull=FALSE,extrapol=T,QP=S)
        W_extra<-W_extra^exposant
        W_extra<-normW(W_extra)
      }
  } else if(method_pred=='model'){
    if(model@Type=='GD') {
      W_extra=kernel_matW(H=c(model@H,model@Ht),kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,diagnull=FALSE,extrapol=T,QP=S)
  } else if(model@Type=='GDT') {
      W_extra=kernel_matW(H=c(model@H,model@Ht),kernels=model@kernels,coords=coords[c(S,O),],NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,diagnull=FALSE,extrapol=T,alpha=model@alpha,QP=S)
      #W_extra=t(apply(W_extra,1,function(x) x*as.numeric(x>=sort(x,decreasing = T)[k_extra+2])))
      W_extra<-W_extra[,tokeep]
      W_extra<-W_extra^exposant
      W_extra<-normW(W_extra)
  } else if(model@Type=='T'){
    W_extra=kernel_matW(H=c(model@Ht),kernels=model@kernels,coords=as.matrix(coords[c(S,O)],ncol=1),NN=model@NN,TP=O,Type=model@Type,adaptive=model@adaptive,diagnull=FALSE,extrapol=T,QP=S)
    }
  } else {
      kernels_pred='shepard'
      adaptive_pred=FALSE
      NNN=k_extra+1
      W_extra=kernel_matW(H=c(k_extra),kernels=kernels_pred,coords=coords,NN=NNN,TP=O,Type=model@Type,adaptive=adaptive_pred,diagnull=FALSE,extrapol=T,QP=S)
    }
}

if(!(model@Model %in% c('OLS','SAR','multiscale_gwr','tds_mgwr','atds_mgwr','atds_gwr')) ) Beta_proj_out=as.matrix(W_extra %*% beta_in[tokeep,])

##############
##### calcul Y_predicted
##############

if(model@Model=='SAR') {
  Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model@X,newdata)),W,e,beta_hat=beta_in,lambda_hat=lambda_in,S,O,type,coords=coords,maxobs=maxobs)

} else if(model@Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv','MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv','MGWRSAR_1_kc_0'))  {
  Y_predicted=BP_pred_SAR(YS,X=as.matrix(rbind(model@X,newdata[,namesX])),W,e,beta_hat=rbind(beta_in,as.matrix(Beta_proj_out)),lambda_hat=c(lambda_in,as.numeric(W_extra %*%lambda_in)),S,O,type,k_extra=k_extra,kernel_extra=kernel_extra,model=model@Model, W_extra = W_extra,coords=coords,maxobs=maxobs)
} else {
  ## ALL VARYING COEF MODEL without spatial autocorrelation
  #browser()
  X=as.matrix(newdata[,namesX])
  Y_predicted=rowSums(as.matrix(Beta_proj_out*X))
}
  }
if(beta_proj) return(list(Y_predicted=Y_predicted,Beta_proj_out=Beta_proj_out))
else return(Y_predicted)
}
