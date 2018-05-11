#' BP_pred_SAR
#' to be documented
#' @usage BP_pred_SAR(YS,X,W,e,beta_hat,lambda_hat,S,O,type='BPN',
#' model='SAR',W_extra=NULL,k_extra=30,kernel_extra='sheppard',coords=NULL)
#' @param YS to be documented
#' @param X to be documented
#' @param W to be documented
#' @param e to be documented
#' @param beta_hat to be documented
#' @param lambda_hat to be documented
#' @param S to be documented
#' @param O to be documented
#' @param type to be documented
#' @param model to be documented
#' @param W_extra to be documented
#' @param k_extra to be documented
#' @param kernel_extra to be documented
#' @param coords to be documented
#' @keywords internal
#' @return to be documented
BP_pred_SAR <-
  function(YS,X,W,e,beta_hat,lambda_hat,S,O,type='BPN',model='SAR',W_extra=NULL,k_extra=30,kernel_extra='sheppard',coords=NULL){
    ###  YS = c(Y_insample)
    ###  X = rbind(X_insample,X_outsample)
    ###  e=residuals_insample
    ###  beta_hat,lambda_hat = coef_insample
    ###  S,O index of insample and outsample

    sigma2=mean(e^2)
    J=S[which(Matrix::rowSums(W[S, O]) > 0)]
    JO=unique(c(sort(J),sort(O)))
    WJ<-W[JO,JO]# WJ<-W[JO,JO]
    rs<-Matrix::rowSums(WJ)
    sup0<-rs>0
    WJ[sup0,]<-WJ[sup0,]/rs[sup0]
    m=length(JO)
    ### creation des index sur J et O
    J2=1:length(J)
    O2=(length(J)+1):m
    if(model!='SAR') {
      if(is.null(W_extra) ) {
        if(ncol(coords)==2) W_extra=KNNX(coords[S,],k_neighbors=k_extra,diagnull=FALSE,kernel=kernel_extra,query=coords[O,]) else {
          stop( "To be coded")
        }
      }
    if(length(lambda_hat)>1) lambda_hat=as.numeric(c(lambda_hat,as.numeric(W_extra %*% as.matrix(lambda_hat))))
    if(length(lambda_hat)>1) lambda_hat<-as.numeric(lambda_hat[JO]) else lambda_hat<-as.numeric(lambda_hat)


    ## Compute Q
    Q=as.matrix((Diagonal(m,1)-lambda_hat*(WJ+Matrix::t(WJ))+lambda_hat^2*(Matrix::t(WJ)%*%WJ))/sigma2)
    ## Compute TC
    if(m<4000) iWJ=solve(as.matrix(Diagonal(m,1)-lambda_hat*WJ)) else if(length(lambda_hat)==1) iWJ=ApproxiW(WJ,lambda_hat,8) else ('Stop predict not possible for varying lambda with more than 4000 obs')
      beta_hat=rbind(beta_hat,as.matrix(W_extra %*% beta_hat))
      YTC = iWJ %*% (as.matrix(rowSums(X[JO,] * beta_hat[JO,])))
    } else {
      Q=as.matrix((Diagonal(m,1)-lambda_hat*(WJ+Matrix::t(WJ))+lambda_hat^2*(Matrix::t(WJ)%*%WJ))/sigma2)
      if(m<4000) iWJ=solve(as.matrix(Diagonal(m,1)-lambda_hat*WJ)) else  iWJ=ApproxiW(WJ,lambda_hat,8)
      YTC= iWJ %*% X[JO,]%*% beta_hat
    }

    ## compute BP ?
    if(type=='BPN') res=YTC[O2]-solve(Q[O2,O2])%*%Q[O2,J2]%*% (YS[J]-YTC[J2])  else res=YTC[O2]
    res
  }
# for(l in 1:max(coords[O,3])) {
#   LindexS<-which(coords[S,3]==l)
#   LindexO<-which(coords[O,3]==l)
#   Wt=KNN(coords[S,][LindexS,],k_extra,diagnull=FALSE,kernel=kernel_extra,query=coords[O,][LindexO,])
#   W_extra[LindexO,LindexS]<-Wt
