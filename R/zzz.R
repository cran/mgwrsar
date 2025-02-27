.onLoad <- function(libname = find.package("mgwrsar"), pkgname = "mgwrsar"){
  # CRAN Note avoidance Thanks to
  if(getRversion() >= "2.15.1")
    utils::globalVariables(c("n","par_model2","regFun2","myseed","prof","spv","ord","myformula_gam","nu","mstop","myblock","Betac_proj_out","Betav_proj_out","Lambdacor","Method","NmaxDist","Penalized","SE","TIME","Type","Wh","Wk","Z","coords","decay","doMC","isgcv","kernels_w","lower_c","lower_cW","lower_d","lower_dW","maxknn","minv","n","n_searchW","ncore","rwild","search_W","upper_c","upper_dW","verbose","xratiomin",'get_ts', 'ALLV', 'Shat', 'TRUEBETA', 'benchmark_GWR', 'get_s', 'i' ,'nns' ,'tol','EvalExpressions' ,'H', 'K' ,'V' ,'alpha' ,'get_Rk', 'model_lm0' ,'n_emp' ,'namesX', 'returned_model', 'theta', 'varying', 'x','plot.mgwrsar'))
  #load(file='data/modc.Rdata')
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("mgwrsar", libpath)
}
utils::globalVariables(c("X", "Y",  "XV", "S", "NN", "MykernelS",  "TP", "k_extra", "dists", "indexG", "Wd", "XC",  "names_betav", "names_betac","ALL_X","KernelTP","tolerance","W", "kernel_w",  "adaptive", "adaptive_W",  "m", "nstage","z","SEV","tS","cell.quadtree","cell.quadtree.leaf","Hp","kWtp",'TP_estim_as_extrapol','eff','kernel_extra','new_XC','new_XV', 'new_data' ,'signif_95','variable'))

if(FALSE){
W_extra=kernel_matW(S=rbind(model@S,newdata_coords),H=mykernels_H,diagnull=FALSE,kernels=mykernels_extra,adaptive=mykernels_adaptive,NN=k_extra[1],Type=model$Type,query_TP=1:nrow(model$S),correctionNB=FALSE,extrapTP=1)
# \item{KernelTP}{A vector containing the kernel types for extrapolation of coefficients from target points. Possible types: shepard ("shepard) default, rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"), gaussian ("gauss"), all in adaptive form.}
# \item{kWtp}{A list of bandwidth size for extrapolation of coefficients from target points.}

S=1:800
O=801:1000
W_extra=kernel_matW(S=coords,H=100,diagnull=FALSE,kernels="shepard",adaptive=FALSE,NN=16,Type='GD',query_TP=S,correctionNB=FALSE,extrapTP=0)
W_extra=W_extra[O,S]
}


