#' assign_control
#' to be documented
#' @usage assign_control(control=list(),n, env = parent.frame())
#' @param control to be documented
#' @param n to be documented
#' @param env env = parent.frame()
#' @keywords internal
#' @return to be documented
assign_control <-
function(control=list(),n, env = parent.frame()) {
#envavirer<<-new.env()
con=list(Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,Method='2SLS',TIME=FALSE,decay=0,minv=1,Type='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,maxknn=500,NmaxDist=5000,Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,remove_local_outlier=FALSE,outv=0,doMC=FALSE,ncore=1,Wh=NULL,SE=FALSE,xratiomin=10e-10)
nmsC <- names(con)
con[(namc <- names(control))] <- control
if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
for(i in 1:length(con))
{
assign(names(con)[i],con[[i]],envir =parent.frame())
#assign(names(con)[i],con[[i]],envir=envavirer)
}
}
