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
    con=list(Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,adaptive=FALSE,Method='2SLS',decay=0,Type='GD',TypeW='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,NN=n,Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,Wh=NULL,SE=FALSE,get_ts=TRUE,get_s=FALSE,get_Rk=FALSE,TP=NULL,KernelTP='Wd',searchB=FALSE,tolerance=0.0001,ncore=1,indexG=NULL,dists=NULL,Wd=NULL,tp_rmse=2,TP_estim_as_extrapol=FALSE,new_data=NULL,new_W=NULL,nu=0.1,mstop=150,family=gaussian(link = "identity"),V=NULL,HRMSE=NULL,HBETA=NULL,alpha=1,criterion='AICc',TP_eval=NULL,kernel_extra=NULL) #

    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
    for(i in 1:length(con))
    {
      assign(names(con)[i],con[[i]],envir =parent.frame())
    }
  }
