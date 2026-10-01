.onLoad <- function(libname, pkgname) {
  .mgwrsar_set_native_threads(1L)
  # CRAN Note avoidance Thanks to
  if(getRversion() >= "2.15.1")
    utils::globalVariables(c("h","X","n","par_model2","regFun2","myseed","prof","spv","ord","myformula_gam","nu","mstop","myblock","Betac_proj_out","Betav_proj_out","Lambdacor","Method","NmaxDist","Penalized","SE","TIME","Type","Wh","Wk","Z","coords","decay","isgcv","kernels_w","lower_c","lower_cW","lower_d","lower_dW","maxknn","minv","n","n_searchW","ncore","rwild","search_W","upper_c","upper_dW","verbose","xratiomin",'get_ts', 'ALLV', 'Shat', 'TRUEBETA', 'benchmark_GWR', 'get_s', 'i' ,'nns' ,'tol','EvalExpressions' ,'H', 'K' ,'V' ,'alpha' ,'get_Rk', 'model_lm0' ,'n_emp' ,'namesX', 'returned_model',  'varying', 'x','plot',"v", "var", "vt", "controlv", "model0", "init_model", "stage1_result","myformula_bk", "data", "coords", "kernels", "Ht", "max_dist", "n", "control","vkst", "k", "i", "init_param_tds", "built_Vseq", "stage1_tds_mgwr", "stage1_tds_mgwr_H", "stage2_atds_mgwr", "stage2_atds_gwr", "%||%",'get_AIC','model','init','maxiter','ratiotime','nstable'))
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("mgwrsar", libpath)
}
utils::globalVariables(c("X", "Y",  "XV", "S", "NN", "MykernelS",  "TP", "k_extra", "dists", "indexG", "Wd", "XC",  "names_betav", "names_betac","ALL_X","KernelTP","tolerance","W", "kernel_w",  "adaptive", "adaptive_W",  "m", "nstage","z","SEV","tS","cell.quadtree","cell.quadtree.leaf","Hp","kWtp",'TP_estim_as_extrapol','eff','kernel_extra','new_XC','new_XV', 'new_data' ,'signif_95','variable'))
