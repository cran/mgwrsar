#' BP_pred_MGWRSARC
#' to be documented
#' @usage BP_pred_MGWRSARC(model,W,S,O,xnewG,varc,coords,type='BPN',k=16)
#' @param model to be documented
#' @param W to be documented
#' @param S to be documented
#' @param O to be documented
#' @param xnewG to be documented
#' @param varc to be documented
#' @param coords to be documented
#' @param type to be documented
#' @param k to be documented
#' @keywords internal
#' @return to be documented
BP_pred_MGWRSARC <-function(model,W,S,O,xnewG,varc,coords,type='BPN',k=16){
    ###  YS = c(Y_insample)
    ###  X = rbind(X_insample,X_outsample)
    ###  e=residuals_insample
    ###  beta_hat,lambda_hat = coef_insample
    ###  S,O index of insample and outsample
    YS = model$Y
    e=model$residuals
    xnew=xnewG[,-which(names(xnew)==varc)]
    X = rbind(model$X,xnewG)

    ### computation of Beta_hat on the whole sample
    wi= as.matrix(GDC((O[1]-1), Y=YS, X=model$X, S=cbind(coords,c(model$S[,3],xnew[,varc])), model$H, model$kernels, minv = 0, maxknn = 500, NmaxDist = 6000,isgcv=FALSE, TIME=FALSE, decay=0))
    for(i in O[-1]) {
      wi= rbind(wi,as.matrix(GDC(i-1, Y=YS, X=model$X, S=cbind(coords,c(model$S[,3],xnew[,varc])), model$H, model$kernels, minv = 0, maxknn = 500, NmaxDist = 6000,isgcv=FALSE, TIME=FALSE, decay=0)))
    }
    Wos=wi[,S]

    ###(Beta_extropolation(coord=coord,k=k,O,type='S+OxS+O'))[-S,-O]
    beta_hat=rbind(beta_hat,as.matrix(Wos %*% beta_hat))
    lambda_hat=model$Betac
    if(length(lambda_hat)>1) lambda_hat=as.numeric(c(lambda_hat,as.numeric(Wos %*% lambda_hat)))
    ###
    sigma2=mean(e^2)
    J=S[which(Matrix::rowSums(W[S, O]) > 0)]
    JO=unique(c(sort(J),sort(O)))
    WJ<-W[JO,JO]
    rs<-Matrix::rowSums(WJ)
    sup0<-rs>0
    WJ[sup0,]<-WJ[sup0,]/rs[sup0]
    m=length(JO)

    if(length(lambda_hat)>1) lambda_hat<-as.numeric(lambda_hat[JO])
    ### creation des index sur J et O
    J2=1:length(J)
    O2=(length(J)+1):m
    ## Compute Q
    Q=(Diagonal(m,1)-lambda_hat*(WJ+Matrix::t(WJ))+lambda_hat^2*(Matrix::t(WJ)%*%WJ))/sigma2
    ## Compute XB
    XB=Matrix::rowSums(X[JO,]* beta_hat[JO,])
    ## Compute TC
    if(m<2000) iWJ=solve(Diagonal(m,1)-lambda_hat*WJ) else iWJ=ApproxiW(WJ,lambda_hat,8)
    YTC=iWJ %*% XB
    ## compute BP
    YBPN=YTC[O2]-solve(Q[O2,O2])%*%Q[O2,J2]%*% (YS[J]-YTC[J2])
    if(type=='BPN') res=YBPN else res=YTC[O2]
    res
  }
