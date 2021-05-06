
#' Computes weight matrix for a given kernel and bandwidth
#'
#' kernelW_C is a function that computes weight matrix for a given kernel and bandwidth
#'
#' @usage kernelW_C(XX, hh, MykernelS, isgcv_, Type, Minv, maxknn_,
#' NmaxDist_, TIME, Decay, DDiagNull)
#' @param XX a matrix with coordinates.
#' @param hh a bandwidth value
#' @param MykernelS a kernel type between ('bin','bisq','gauss','gauss_adapt',
#' 'knn',bisq_knn')
#'@param isgcv_ default FALSE for computing CV criteria (for example for selecting optimal bandwidth).
#' @param Type  Kernel type.
#' @param Minv, minimum number of neighbors when using distance kernels.
#' @param maxknn_ default 500, when n>NmaxDist only  maxknn first neighbours are used for computation distance
#' @param NmaxDist_ default 5000, when n>NmaxDist only  maxknn first neighbours are used for computing distance
#' @param TIME default FALSE, time is used for computing weigths if TIME is TRUE
#' weigth for future observation are set to zero
#' @param Decay time decay when time is used for computing weigths.
#' @param DDiagNull default FALSE, if TRUE diagonal has zero weights.
#' @return  kernelW_C returns a sparse weight matrix
#' @seealso  MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar
#' @examples
#' \donttest{
#' data(mydata);coord=as.matrix(mydata[,c("x_lat","y_lon")]);
#' W=kernelW_C(coord,100,'bisq_knn',FALSE,'GD',1,500,5000,FALSE,0,FALSE)
#' plot(D_dense_C(coord[1,1],coord[1,2],coord[,1],coord[,2]),W[1,])
#' }
kernelW_C <-
function (XX, hh, MykernelS, isgcv_, Type, Minv, maxknn_, NmaxDist_,
    TIME, Decay, DDiagNull)
.Call("kernelW_C", XX, hh, MykernelS, isgcv_, Type, Minv, maxknn_,
    NmaxDist_, TIME, Decay, DDiagNull, PACKAGE = "mgwrsar")
