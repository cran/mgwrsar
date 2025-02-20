#' dvcov
#' return diagonal of vcov of lm.fit object
#' @usage dvcov(lmfit)
#' @param lmfit  a lm.fit object
#' @noRd
#' @return to be documented
dvcov = function(lmfit) {
  qr=lmfit$qr
  residuals=lmfit$residuals
  p <- qr$rank
  n <-length(residuals)
  qr.mat <- as.matrix(qr$qr[1L:p, 1L:p])
  qr.mat[row(qr.mat) > col(qr.mat)] <- 0
  qrinv <- solve(qr.mat)
  s2 = sum(residuals^2) / (n-p)
  vcov <- qrinv %*% t(qrinv) * s2
  diag(vcov)
}
