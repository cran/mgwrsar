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
    con=list(Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,adaptive=F,Method='2SLS',TIME=FALSE,decay=0,Type='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,NN=min(2000,n),Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,remove_local_outlier=FALSE,outv=0,Wh=NULL,SE=FALSE,TP=NULL,kWtp=16,KernelTP='Wd',nstop=NULL,nneg=8,Wtp=NULL,searchB=F,tolerance=0.0001,doMC=FALSE,ncore=1,indexG=NULL,dists=NULL,Wd=NULL,eta=1,tp_rmse=2,S_out=NULL,new_data=NULL,new_W=NULL,pred=NULL) #

    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
    for(i in 1:length(con))
    {
      assign(names(con)[i],con[[i]],envir =parent.frame())
    }
  }
