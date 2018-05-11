.onLoad <- function(libname = find.package("mgwrsar"), pkgname = "mgwrsar"){
  # CRAN Note avoidance Thanks to
  if(getRversion() >= "2.15.1")
    utils::globalVariables(c("Lambdacor","Method","NmaxDist","Penalized","SE","TIME","Type","Wh","Wk","Z","coord","decay","doMC","isgcv","kernels_w","lower_c","lower_cW","lower_d","lower_dW","maxknn","minv","n","n_searchW","ncore","outv","remove_local_outlier","rwild","search_W","upper_c","upper_dW","verbose","xratiomin")) #,"envavirer"
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("mgwrsar", libpath)
}
