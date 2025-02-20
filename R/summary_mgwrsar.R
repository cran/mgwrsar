#' Print a summary of mgwrsar models
#' @usage summary_mgwrsar(model)
#' @param model a model of class mgwrsar
#' @return a summary of mgwrsar models
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 4 nearest neighbors with 0 in diagonal
#'  W=kernel_matW(H=4,kernels='rectangle',coords=coords,NN=4,adaptive=TRUE,
#'  diagnull=TRUE)
#'  mgwrsar_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,
#'  coords=coords, fixed_vars='X2',kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_kc_kv',
#'  control=list(SE=FALSE,adaptive=TRUE,W=W))
#'  summary_mgwrsar(mgwrsar_0_kc_kv)
#' }
#' @noRd
summary_mgwrsar <- function(model) {
  if(!is(model,'mgwrsar')) stop("not a mgwrsar object")
  cat("Call:\n")
  print(model$mycall)
  n <- length(model$Y)
  cat("Model:", model$Model, "\n")
  if(!(model$Model %in% c('OLS','GWR','MGWR'))) cat("Method for spatial autocorrelation:", model$Method, "\n")
  cat("Kernels function:", model$kernels, "\n")
  cat("Kernels adaptive:", ifelse(model$adaptive,'YES','NO'), "\n")
  cat("Kernels type:", model$Type, "\n")
  cat("Bandwidth:", model$H, "\n")
  cat("Computation time:", model$ctime, "\n")
  cat("Use of parallel computing:", model$doMC, " ncore =",model$ncore,"\n")
  cat("Use of rough kernel:",ifelse(model$NN<n & (!model$adaptive | (model$adaptive & model$kernels=='gauss')) ,paste0('YES, ',model$NN,' neighbors / ',n),'NO'),"\n")
  cat("Use of Target Points:", ifelse(is.null(model$TP) | length(model$TP)==n,'NO','YES'), "\n")
   if(length(model$TP)!=n) {
  # cat(paste0("Use of ",ifelse(length(model$TP)>1,length(model$K),'single ')," pass of Target Points \n"))
  # cat("Use of iterative mode of TP choice : ",ifelse(model$TPboost,'yes','no')," \n")
   cat("Number of Target Points ", length(unlist(model$TP)), "\n")
   }
  cat("Number of data points:", n, "\n")
  if(!is.null(model$XC)) {cat("Constant parameters:", names(model$Betac), "\n")
  cat(as.numeric(model$Betac),'\n')}
  if(!is.null(model$XV)){ cat("Varying parameters:", colnames(model$XV), "\n")
  CM <- apply(model$Betav, 2, summary)
  printCoefmat(CM)
  }
  if (length( model$edf)>0) {
    cat("Effective degrees of freedom:", model$edf, "\n")
  }
  if (!is.null( model$tS)) cat("AICc:", model$AICc, "\n")
  if (!is.null(model$CV)) cat("LOOCV:", model$CV, "\n")

  cat("Residual sum of squares:", model$SSR, "\n")
  cat("RMSE:", model$RMSE, "\n")

  invisible(model)
}
