#' Title GDXC5
#' to be documented
#' @usage GDXC5(j, Y, X, S, H, kernels, minv = 0, maxknn = 500, NmaxDist = 6000,isgcv, TIME, decay)
#' @param j   to be documented
#' @param Y  to be documented
#' @param X  to be documented
#' @param S  to be documented
#' @param H  to be documented
#' @param kernels  to be documented
#' @param minv  to be documented
#' @param maxknn  to be documented
#' @param NmaxDist  to be documented
#' @param isgcv  to be documented
#' @param TIME  to be documented
#' @param decay  to be documented
#' @keywords internal
#' @return to be documented
GDXC5 <-
function(j, Y,X,S, H, kernels, minv = 0, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay){
  mat<-matrix(1,5,5)
  mat[1,2]<-H[3]
  mat[1,3]<-H[4]
  mat[1,4]<-H[5]
  mat[1,5]<-H[6]
  mat[2,1]<-H[7]
  mat[2,3]<-H[8]
  mat[2,4]<-H[9]
  mat[2,5]<-H[10]
  mat[3,1]<-H[11]
  mat[3,2]<-H[12]
  mat[3,4]<-H[13]
  mat[3,5]<-H[14]
  mat[4,1]<-H[15]
  mat[4,2]<-H[16]
  mat[4,3]<-H[17]
  mat[4,5]<-H[18]
  mat[5,1]<-H[19]
  mat[5,2]<-H[20]
  mat[5,3]<-H[21]
  mat[5,4]<-H[22]
  Wz<-mat[S[j + 1, 4],S[, 4]]
 kernel_C(as.matrix(S[, 1:2]), j, H[1], kernels[1], minv, FALSE, decay, isgcv, TRUE)*kernel_C(as.matrix(S[, 3]), j, H[2], kernels[2], minv, FALSE, decay, isgcv, TRUE)*Wz
}
