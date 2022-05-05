#' BP_pred_SAR
#' to be documented
#' @usage BP_pred_SAR(YS,X,W,e,beta_hat,lambda_hat,S,O,type='BPN',
#' model='SAR',W_extra=NULL,k_extra=30,kernel_extra='sheppard',coords=NULL,maxobs=4000)
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
#' @param maxobs to be documented
#' @noRd
#' @return to be documented
BP_pred_SAR<-function (YS, X, W, e, beta_hat, lambda_hat, S, O, type = "BPN",model = "SAR", W_extra = NULL, k_extra = 30, kernel_extra = "sheppard", coords = NULL,maxobs=4000){
  warnings('YS  and e should have S rows, X should be the matrix of new_data with S+O rows (order S then O), W should be a row normalized matrix of dim (S+0,S+O)')
  n=length(c(S,O))
  if(n!=nrow(X)) stop('Wrong arguments: S should be the row index of in-sample and O the row index of out-sample, ex. S=1:100,O=101:120 ')
  sigma2 = mean(e^2)
  if(n < maxobs) iW = solve(as.matrix(Diagonal(n, 1) - lambda_hat *W)) else iW = ApproxiW(W, lambda_hat, 8)
  if(is.matrix(beta_hat))  YTC = iW %*% rowSums(X* beta_hat)  else  YTC = iW %*% X %*% beta_hat
  if (type == "BPN"){
    Q = (Diagonal(n, 1) - lambda_hat * (W + Matrix::t(W)) + lambda_hat^2 * (Matrix::t(W) %*% W))/sigma2
    res = YTC[O] - solve(Q[O, O]) %*% Q[O, S] %*% (YS[S] -YTC[S])
  }
else res = YTC[O]
  as.numeric(res)
}
