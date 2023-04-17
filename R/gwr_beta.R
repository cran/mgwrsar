#' gwr_beta
#' to be documented
#' @usage gwr_beta(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=F,SE=FALSE,
#' remove_local_outlier=F,outv=0.01,KernelTP='sheppard',kWtp=8,
#' doMC=FALSE,ncore=1)
#' @param Y A vector of response
#' @param XV A matrix with covariates with non stationnary parameters
#' @param ALL_X A matrix with all covariates
#' @param coord A matrix with spatial coordinates
#' @param TP An index of target points.
#' @param indexG Precomputed Matrix of indexes of NN neighbors.
#' @param Wd Precomputed Matrix of weights.
#' @param NN Number of spatial Neighbours for kernels computations
#' @param W The spatial weight matrix for spatial dependence
#' @param isgcv leave one out cross validation, default FALSE
#' @param SE If standard error are computed, default FALSE
#' @param remove_local_outlier  Remove local outlier
#' @param outv percentile threshold for outlier.
#' @param KernelTP  Kernel type for extrapolation of Beta from Beta(TP)
#' @param kWtp  Number of neighbours for extrapolation of Beta from Beta(TP)
#' @param doMC  Boolean for parallel computation.
#' @param ncore  Number of cores for parallel computation.
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=F,SE=FALSE,remove_local_outlier=F,outv=0.01,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,pred=FALSE)
{
  n=length(Y)
  ntp=length(TP)
  if(pred) {
    SE=FALSE
    isgcv=F
  }
  if(!is.null(XV)) m=ncol(XV) else m=0
  tS=0
  namesXV=colnames(XV)
  if (!is.null(W)) {
    PhWy=PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,n))
    XV = cbind(XV,PhWy)
  }## ce point pose probleme dans le cas isgcv

  if(isgcv) loo=-1 else loo=1:NN
  if(doMC) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP))
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    if(SE) {
      tS<-0;
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
    }
    if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
    for(z in myblocks[[myblock]]){
      index=indexG[z,loo] #### commencer ici les adaptations boost
      wd<-sqrt(Wd[z,loo])
      Yw<-wd*Y[index]
      Xw=wd*XV[index,]
      if(remove_local_outlier){
        cooks=cooks.distance(lm(Yw~Xw-1))
        i_reject=which(cooks>quantile(cooks,1-outv))
        Xw<-Xw[-i_reject,]
        Yw<-Yw[-i_reject]
      }
      lml=lm.fit(as.matrix(Xw),as.matrix(Yw))
      betav=lml$coefficients ### finir ici pour glmboost
      coefNA<-which(is.na(betav))
      betav[coefNA]<-0
      if(length(coefNA)>0) {
        lml=lm.fit(as.matrix(Xw[,-coefNA]),as.matrix(Yw))
        betav[-coefNA]<-lml$coefficients
      }  ### finir ici pour glmboost
      if(!is.null(W) & abs(betav[m + 1])>1) {
        betav[m + 1]=sign(betav[m + 1])*0.99
        if(m>0) {
          lml=lm.fit(as.matrix(Xw[,-c(coefNA,m+1)]),as.matrix(Yw-betav[m + 1]*Xw[,m + 1]))
          betav[setdiff(1:m,coefNA)]=lml$coefficients
        }
      }
      if(SE & !isgcv) {
        coef_NON_NA=setdiff(1:ncol(Xw),coefNA)
        rss <- sum(lml$residuals^2)
        rdf <- length(Yw) - ncol(Xw)+length(coefNA)
        resvar <- rss/rdf
        R <- chol2inv(lml$qr$qr)
        diagR=diag(R)
        SEV[TP[z],coef_NON_NA] <- sqrt(diagR * resvar)
        Zwi = try(solve(crossprod(Xw[, coef_NON_NA], Xw[, coef_NON_NA])) %*% t(Xw[, coef_NON_NA]),silent = TRUE)
        tS=tS+ifelse(class(Zwi)[1]=='try-error',0,(Xw[1, coef_NON_NA] %*% Zwi)[, 1])
      }
      if(!pred){ Betav[TP[z],]<-betav } else {Betav[z,]<-betav}
    }
    if(!(SE & !isgcv)) {
      sev=NULL
      tS=NULL
    } else {
      sev=SEV[TP[myblocks[[myblock]]],]
    }
    rm(index,wd,Yw,Xw,betav)
    gc()
    list(betav=Betav[TP[myblocks[[myblock]]],],sev=sev,tS=tS)
  } #

  if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))

  if(!pred) {
    Betav[TP,]<-res$betav
    if(SE) {
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
      SEV[TP,]=res$sev
      tS=sum(res$tS)
    }
  } else {
    Betav<-res$betav
  }

  if(SE) colnames(SEV)=colnames(XV)
  if(ntp<length(Y) & !pred){
    #if(KernelTP=='Wd'){ ### on ne recalcule pas W
    Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
    #}

    # else if(KernelTP=='W1') {
    #   ### on recalcule W same kernel/bandwdith post - normalisation
    #   Wtp=kernel_matW(S=coord,H=H,diagnull=FALSE,kernels=kernels,adaptive=adaptive,NN=NN,Type='GD',query_TP=NULL,rowNorm=FALSE,correctionNB=FALSE,extrapTP=0)[-TP,TP]
    #  Wtp=Wtp/Matrix::rowSums(Wtp)
    #  } else {  ### On utilise un autre kernel/bandwidth
    #   Wtp=kernel_matW(S=coord,H=kWtp,diagnull=FALSE,kernels='sheppard',adaptive=FALSE,NN=kWtp,Type='GD',query_TP=TP,rowNorm=TRUE,correctionNB=FALSE,extrapTP=1)[-TP,TP]
    # }
    Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
    if(SE) SEV[-TP,]=as.matrix(Wtp%*% SEV[TP,])
  }
  if(is.null(W))  colnames(Betav)=namesXV else colnames(Betav)=c(namesXV,'lambda')
  if(SE & !isgcv & !pred) list(Betav=Betav,SEV=SEV,edf=n-tS,tS=tS) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL)
}



# gwr_beta_derecated<-function (Y, XV, ALL_X, TP, indexG, Wd, NN, W = NULL, isgcv = F,
#           SE = FALSE, remove_local_outlier = F, outv = 0.01, kernels = NULL,
#           H = NULL, adaptive = NULL, doMC = FALSE, ncore = 1, pred = FALSE)
# {
#   n = length(Y)
#   ntp = length(TP)
#   if (pred) {
#     SE = FALSE
#     isgcv = F
#   }
#   if (!is.null(XV))
#     m = ncol(XV)
#   else m = 0
#   if (SE) {
#     tS <- 0
#     edf = n
#     SEV <- matrix(0, nrow = n, ncol = ifelse(is.null(W),
#                                              m, m + 1))
#   }
#   namesXV = colnames(XV)
#   if (!is.null(W)) {
#     PhWy = PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,
#                                                          n))
#     XV = cbind(XV, PhWy)
#   }
#   if (pred)
#     Betav = matrix(0, nrow = ntp, ncol = ifelse(is.null(W),
#                                                 m, m + 1))
#   else Betav = matrix(0, nrow = n, ncol = ifelse(is.null(W),
#                                                  m, m + 1))
#   if (isgcv)
#     loo = -1
#   else loo = 1:NN
#   if (doMC) {
#     registerDoParallel(cores = ncore)
#   }
#   else registerDoSEQ()
#   for (z in 1:length(TP)) {
#     index = indexG[z, loo]
#     wd <- sqrt(Wd[z, loo])
#     Yw <- wd * Y[index]
#     Xw = wd * XV[index, ]
#     if (remove_local_outlier) {
#       cooks = cooks.distance(lm(Yw ~ Xw - 1))
#       i_reject = which(cooks > quantile(cooks, 1 - outv))
#       Xw <- Xw[-i_reject, ]
#       Yw <- Yw[-i_reject]
#     }
#     lml = lm.fit(as.matrix(Xw), as.matrix(Yw))
#     betav = lml$coefficients
#     coefNA <- which(is.na(betav))
#     betav[coefNA] <- 0
#     if (length(coefNA) > 0) {
#       lml = lm.fit(as.matrix(Xw[, -coefNA]), as.matrix(Yw))
#       betav[-coefNA] <- lml$coefficients
#     }
#     if (!is.null(W) & abs(betav[m + 1]) > 1) {
#       betav[m + 1] = sign(betav[m + 1]) * 0.99
#       if (m > 0) {
#         lml = lm.fit(as.matrix(Xw[, -c(coefNA, m + 1)]),
#                      as.matrix(Yw - betav[m + 1] * Xw[, m + 1]))
#         betav[setdiff(1:m, coefNA)] = lml$coefficients
#       }
#     }
#     if (SE & !isgcv) {
#       coef_NON_NA = setdiff(1:ncol(Xw), coefNA)
#       rss <- sum(lml$residuals^2)
#       rdf <- length(Yw) - ncol(Xw) + length(coefNA)
#       resvar <- rss/rdf
#       R <- chol2inv(lml$qr$qr)
#       diagR = diag(R)
#       SEV[TP[z], coef_NON_NA] <- sqrt(diagR * resvar)
#       Zwi = try(solve(crossprod(Xw[, coef_NON_NA], Xw[,
#                                                       coef_NON_NA])) %*% t(Xw[, coef_NON_NA]), silent = TRUE)
#       tS = tS + ifelse(class(Zwi)[1] == "try-error", 0,
#                        (Xw[1, coef_NON_NA] %*% Zwi)[, 1])
#     }
#     else {
#       sev = NULL
#       tS = NULL
#     }
#     if (!pred) {
#       Betav[TP[z], ] <- betav
#     }
#     else {
#       Betav[z, ] <- betav
#     }
#   }
#   if (SE)
#     colnames(SEV) = colnames(XV)
#   if (ntp < length(Y) & !pred) {
#     Wtp <- normW(Matrix::t(sparseMatrix(i = rep(1:ntp, each = NN),
#                                         j = as.numeric(t(indexG)), dims = c(ntp, n), x = as.numeric(t(Wd))))[-TP,
#                                         ])
#     Betav[-TP, ] = as.matrix(Wtp %*% Betav[TP, ])
#     if (SE)
#       SEV[-TP, ] = as.matrix(Wtp %*% SEV[TP, ])
#   }
#   if (is.null(W))
#     colnames(Betav) = namesXV
#   else colnames(Betav) = c(namesXV, "lambda")
#   if (SE & !isgcv & !pred)
#     list(Betav = Betav, SEV = SEV, edf = n - tS, tS = tS)
#   else list(Betav = Betav, SEV = NULL, edf = NULL, tS = NULL)
# }
