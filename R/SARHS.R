#' SARHS
#' to be documented
#' @usage SARHS(Xw,YYw,SE,W,Betav)
#' @param Xw  to be documented
#' @param YYw  to be documented
#' @param SE  to be documented
#' @param W  to be documented
#' @param Betav  to be documented
#' @keywords internal
#' @return to be documented
SARHS <-
function(Xw,YYw,SE,W,Betav){
m<-length(Betav)
rho=sign(Betav[m])*.99
YYw<-YYw-rho*W%*%YYw
modelsar<-fastlmLLT_C(as.matrix(Xw[,-m]),as.matrix(YYw),SE)
list(Betav=c(modelsar$Betav,rho),se=modelsar$se)
}
