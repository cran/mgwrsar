#' assign_control
#' to be documented
#' @usage assign_control(control=list(),n, env = parent.frame())
#' @param control to be documented
#' @param n to be documented
#' @param env env = parent.frame()
#' @noRd
#' @return to be documented
assign_control <-
  function(control=list(),n, env = parent.frame()) {
    con=list(Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,adaptive=FALSE,Method='2SLS',TIME=FALSE,decay=0,Type='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,NN=min(2000,n),Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,Wh=NULL,SE=FALSE,get_ts=FALSE,TP=NULL,kWtp=16,KernelTP='Wd',nstop=NULL,nneg=8,Wtp=NULL,searchB=FALSE,tolerance=0.0001,doMC=FALSE,ncore=1,indexG=NULL,dists=NULL,Wd=NULL,eta=1,tp_rmse=2,S_out=FALSE,new_data=NULL,new_W=NULL,pred=NULL,nu=0.1,mstop=150,family=gaussian(link = "identity")) #

    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
    for(i in 1:length(con))
    {
      assign(names(con)[i],con[[i]],envir =parent.frame())
    }
  }
