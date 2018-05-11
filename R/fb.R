#' fb
#' to be documented
#' @usage fb(formula,data,coord,fixed_vars,Model,control,kernel,e_search)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param kernel to be documented
#' @param e_search to be documented
#' @keywords internal
#' @return to be documented
fb <-
function(formula,data,coord,fixed_vars,Model,control,kernel,e_search){
			control$isgcv<-TRUE
			starting=10^100
			if(e_search$search_W){
			cat('\n########\nSearch W stage\n########')
			for(kernel_w in e_search$kernels_w){
			cat(paste('\n########\nKernel candidate for W : ',kernel_w,'\n########',sep=''))
			if(e_search$search_W & !(Model %in% c('GWR','MGWR'))) {
				control$W<-NULL
                if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv','MGWRSAR_1_kc_0')) {
                	e_search=support_search(formula, data,  coord, fixed_vars=NULL, Model='GWR',control,kernel,e_search)
                	opt=bandwidth_search(formula, data,  coord, fixed_vars=NULL, Model='GWR',control,kernel,e_search)
                } else if(Model %in% c('MGWRSAR_1_kc_kv','MGWRSAR_0_kc_kv')){
               		e_search=support_search(formula, data,  coord, fixed_vars, Model='MGWR',control,kernel,e_search)
                	opt=bandwidth_search(formula, data,  coord, fixed_vars, Model='MGWR',control,kernel,e_search)
                }
                if(Model!='SAR') e_search$h=opt$minimum else e_search$h=0
                if(kernel_w=='bisq' & Model!='SAR'){
               control2=control
               control2$W<-kernelW_C(coord,e_search$lower_cW,kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE)  #control$Type --> Type
               while(is.na(SSR_h(formula, data,  coord, fixed_vars,kernel, H=e_search$h, Model, control=control2,e_search$Penalized)) & e_search$lower_cW<e_search$upper_cw){
               e_search$lower_cW=e_search$lower_cW*1.2
               control2$W<-kernelW_C(coord,e_search$lower_cW,kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
               }
               }
               optW=bandwidth_searchW(formula, data,  coord,fixed_vars,Model,control,kernel,kernel_w,e_search)
				control$W<-kernelW_C(coord,optW$minimum,kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
				if(n_searchW>1 & Model!='SAR'){
				for(i in 2:n_searchW){
				cat(paste('\n## Round ',i,'/',e_search$n_searchW,'\n',sep=''))
				cat('optW$objective = ',optW$objective,'optW$minimum = ',optW$minimum,' opt$objective = ',opt$objective,' opt$minimum = ',opt$minimum,'\n')
				opt=bandwidth_search(formula, data,  coord, fixed_vars, Model=Model,control,kernel,e_search)
				e_search$h=opt$minimum
				optW=bandwidth_searchW(formula, data,  coord,fixed_vars,Model,control,kernel,kernel_w,e_search)
				control$W<-kernelW_C(coord,optW$minimum,kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
				}
				}
				}
				if(optW$objective<starting) {
 				e_search$kernel_w=kernel_w ### new
 				e_search$h_w=optW$minimum ### new
				starting=optW$objective}
				}
			control$W<-kernelW_C(coord,e_search$h_w,e_search$kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
			}
		if(Model!='SAR'){
		cat('\n########\nSearch bandwidth stage\n########\n')
		e_search=support_search(formula, data,  coord, fixed_vars, Model,control,kernel,e_search)
        opt=bandwidth_search(formula, data,  coord, fixed_vars, Model,control,kernel,e_search)
        cat(' opt$objective = ',opt$objective,' opt$minimum = ',opt$minimum,'\n')
        control$isgcv=FALSE
        config_model<-c(Model,kernel[1],opt$minimum,e_search$kernel_w,e_search$h_w)
        } else config_model<-c(Model,'unknown',0,e_search$kernel_w,e_search$h_w)
 fill_DGPTAB(formula,data,coord,fixed_vars,Model,control,opt,e_search$search_W,config_model)
}
