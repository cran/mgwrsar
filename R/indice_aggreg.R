#' indice_aggreg
#' to be documented
#' @usage indice_aggreg(D)
#' @param D  to be documented
#' @noRd
#'
#' @return to be documented
indice_aggreg <-
function(D){
#rx=max(coords[,1])-min(coords[,1])
#ry=max(coords[,2])-min(coords[,2])
#surf=rx*ry
surf=1
lambda=surf/n
dmin=apply(D,1,function(x) sort(x)[2])
ed=sqrt(lambda)
dbar=mean(abs(dmin-ed))
(dbar/ed-0.5)*2  ### verifier si 0.5 est independant de S et n
}
