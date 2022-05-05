#' Sidak_cor_MGWR
#' to be documented
#' @usage Sidak_cor_MGWR(alpha=0.05,modelH0,modelH1,D=NULL,type='rn')
#' @param alpha to be documented
#' @param modelH0 to be documented
#' @param modelH1 to be documented
#' @param D to be documented
#' @param type to be documented
#' @noRd
#' @return to be documented
Sidak_cor_MGWR <-
function(alpha=0.05,modelH0,modelH1,D=NULL,type='rn') {
if(is.null(D) & type=='df0/edf1*r*2*n*(1-ia)^4') stop("When type='df0/edf1*r*2*n*(1-ia)^4' you have to provide D a nxn distance matrix")
n<-nrow(modelH0$residuals)
if(type=='rn') r=(ncol(modelH1$XV)-ncol(modelH0$XV))*n
else if(type=='n') r=n else if(type=='edf0') r=modelH0$edf
else if(type=='edf1') r=modelH1$edf
else if(type=='r*edf0') r=(ncol(modelH1$XV)-ncol(modelH0$XV))*modelH0$edf
else if(type=='n+(edf0-edf1)') r=(n+(modelH0$edf-modelH1$edf))
else if(type=='edf0-edf1') r=(modelH0$edf-modelH1$edf)
else if(type=='edf0/edf1*r*2*n') r=2*modelH0$edf/modelH1$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))*n
else if(type=='edf0/edf1*r*2*n*(1-ia)^4') {
ia<-indice_aggreg(D)
r=2*modelH0$edf/modelH1$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))*n*(1-ia)^4
}
else if(type=='edf0/edf1*r*2*n*(1-ia)^3') {
ia<-indice_aggreg(D)
r=2*modelH0$edf/modelH1$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))*n*(1-ia)^3
}
else if(type=='edf0/edf1*r*2*n*(1-ia)^2') {
ia<-indice_aggreg(D)
r=2*modelH0$edf/modelH1$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))*n*(1-ia)^2
}
else if(type=='2*(edf0-edf1)*r') {
r=2*(modelH0$edf-modelH1$edf)*(ncol(modelH1$XV)-ncol(modelH0$XV))
}
else if(type=='(edf0-edf1)*r') {
r=2*(modelH0$edf-modelH1$edf)*(ncol(modelH1$XV)-ncol(modelH0$XV))
}
else if(type=='r*edf0') {
r=modelH0$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))
}
else if(type=='r*edf1') {
r=modelH1$edf*(ncol(modelH1$XV)-ncol(modelH0$XV))
}
1-(1-alpha)^(1/(r))
}
