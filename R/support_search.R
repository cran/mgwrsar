#' support_search
#' to be documented
#' @usage support_search(formula, data,  coord, fixed_vars, Model,control,kernel,e_search)
#' @param formula  to be documented
#' @param data  to be documented
#' @param coord  to be documented
#' @param fixed_vars  to be documented
#' @param Model  to be documented
#' @param control  to be documented
#' @param kernel  to be documented
#' @param e_search  to be documented
#' @keywords internal
#' @return to be documented
support_search<-function(formula, data,  coord, fixed_vars, Model,control,kernel,e_search){
  cat('Bandwidth h: Searching support for kernel',kernel,' with Model',Model,'\n');
  if(kernel[1] %in% c('bisq','trisq','gauss','bin','epanechnikov','bisq_ext','gauss_ext')) {
    ### lower
    starting<-e_search$lower_cw

    change=TRUE
    while(change) {
      temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$lower_cw, Model, control=control,e_search$Penalized),silent=TRUE)
      change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$lower_cw<e_search$upper_cw
      if(change) e_search$lower_cw=e_search$lower_cw*1.2;}


    if(e_search$lower_cw>starting) {
      e_search$lower_cw=e_search$lower_cw/1.2;
      change=TRUE
      while(change) {
        temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$lower_cw, Model, control=control,e_search$Penalized),silent=TRUE)
        change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$lower_cw<e_search$upper_cw
        if(change) e_search$lower_cw=e_search$lower_cw*1.01;
        }
    }
    if(e_search$lower_cw>starting) cat('Increase lower bound lower_c to',e_search$lower_cw,' \n');

    #### upper
    starting<-e_search$upper_cw
    change=TRUE
    while(change) {
      temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$upper_cw, Model, control=control,e_search$Penalized),silent=TRUE)
      change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$upper_cw>e_search$lower_cw
      if(change) e_search$upper_cw=e_search$upper_cw*0.8;
      }

    if(e_search$upper_cw<starting) {
      e_search$upper_cw=e_search$upper_cw/0.8;
      change=TRUE
      while(change) {
        temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$upper_cw, Model, control=control,e_search$Penalized),silent=TRUE)
        change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$upper_cw>e_search$lower_cw
        if(change) e_search$upper_cw=e_search$upper_cw*0.99;}
    }
    if(e_search$upper_cw<starting) cat('Decrease upper bound upper_c to',e_search$upper_cw,' \n');
  } else if(kernel[1] %in% c('bisq_knn','gauss_adapt','gauss_knn','knn','epanechnikov_knn')) {

    ### lowerd
    starting<-e_search$lower_dw
    change=TRUE
    while(change) {
      temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$lower_dw, Model, control=control,e_search$Penalized),silent=TRUE)
      change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$lower_dw<0.65*e_search$n
      if(change) e_search$lower_dw=e_search$lower_dw*1.05;
      }

    if(e_search$lower_dw>starting) {
      e_search$lower_dw=e_search$lower_dw/1.05;
      change=TRUE
      while(change) {
        temp_ssr=try(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$lower_dw, Model, control=control,e_search$Penalized),silent=TRUE)
        change=(is.na(temp_ssr)|class(temp_ssr)=='try-error') & e_search$lower_dw<0.65*e_search$n
        if(change) e_search$lower_dw=ceiling(e_search$lower_dw*1.01);
        }
    }
    if(e_search$lower_dw>starting) cat('Increase lower bound lower_d to',e_search$lower_dw,' \n');} else {stop('Unknown kernel\n')}
  e_search
}

