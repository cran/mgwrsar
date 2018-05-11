#' Print a summary of mgwrsar models
#' @usage summary_mgwrsar(model)
#' @param model a model of class mgwrsar
#' @return a summary of mgwrsar models
#' @seealso  MGWRSAR, bandwidths_mgwrsar, plot_mgwrsar, predict_mgwrsar, kernelW_C
#'
#' @examples
#' \dontrun{
#' library(mgwrsar)
#' data(data_mgwrsar)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' W=KNN(coord,8)
#' model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#'  fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model = 'GWR',
#'  control=list(SE=TRUE))
#' summary_mgwrsar(model_GWR)
#' }
summary_mgwrsar <- function(model) {
  if(class(model) != "mgwrsar") stop("not a mgwrsar object")
  cat("Call:\n")
  print(model$mycall)
  n <- length(model$Y)
  cat("Model:", model$Model, "\n")
  if(!(model$Model %in% c('OLS','GWR','MGWR'))) cat("Method for spatial autocorrelation:", model$Model, "\n")
  cat("Kernels function:", model$kernels, "\n")
  cat("Kernels type:", model$type, "\n")
  cat("bandwidth:", model$H, "\n")
  cat("Number of data points:", n, "\n")
  if(!is.null(model$XC)) {cat("Constant parameters:", colnames(model$XC), "\n")
  cat(as.numeric(model$Betac),'\n')}
  if(!is.null(model$XV)){ cat("Varying parameters:", colnames(model$XV), "\n")
  CM <- apply(model$Betav, 2, summary)
  printCoefmat(CM)
  }
  if (length( model$edf)>0) {
    cat("Effective degrees of freedom:", model$edf, "\n")
    cat("AIC", 2*n*log(sqrt(model$SSR/n)) + n*log(2*pi) + n + model$tS, "\n")
  }
  cat("Residual sum of squares:", model$SSR, "\n")
  invisible(model)
}
