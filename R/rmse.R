#' Root Mean Square Error
#'
#' Computes the RMSE of a vector: \code{sqrt(mean(err^2))}.
#'
#' @usage rmse(err)
#' @param err Numeric vector of errors or residuals.
#' @noRd
#' @return A scalar RMSE value.
rmse<-function(err) sqrt(mean(err^2))
