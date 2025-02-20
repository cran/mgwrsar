#' Asymptotic test for mgwrsar models using Sidak correction,
#' EXPERIMENTAL (alpha=0.05 percent)
#' @usage Tf(model_A,model_B,pv=FALSE)
#' @param model_A  A mgwrsar model H0
#' @param model_B  A mgwrsar model H1
#' @param pv  If TRUE p-value are estmimated using Sidak correction, default FALSE
#' @noRd
#' @return to be documented
Tf <-
function(model_A,model_B,pv=FALSE){
# if(model_A$Model=='MGWRSAR_0_0_kv') {
# 	model_A=model_B
# 	model_B=model_A
# } else if (model_A$Model=='GWR') {
# 	if(model_B$Model=='MGWRSAR_0_0_kv'){
# 		model_A=model_A
# 		model_B=model_B
# 	}
# 	else {
# 		model_A=model_B
# 		model_B=model_A
# 	}
# } else if (model_A$Model=='SAR') {
# 	if(model_B$Model=='MGWRSAR_0_0_kv'){
# 		model_A=model_A
# 		model_B=model_B
# 		}
# 	else if(model_B$Model=='OLS') {
# 		model_A=model_B
# 		model_B=model_A
# 	} else if(model_B$Model=='GWR') stop('NA')
# } else if (model_A$Model=='OLS') {
# model_A=model_A
# model_B=model_B
# }
  n=length(model_A@Y)
T<-( (model_A@SSR-model_B@SSR)/(model_A@edf-model_B@edf) )/ (model_B@SSR/model_B@edf)

    if(pv==TRUE) {
        if(model_A@Model %in% c('MGWR') & model_B@Model=='GWR') {
        	pv<-1-pf(T, (model_B@tSS-model_A@tSS)^2/(model_B@tSS2-model_A@tSS2), (n-model_B@tSS)^2/(n-model_B@tSS2))
        	Accept_HO=(pv<Sidak_cor_MGWR(alpha=0.05,model_A,model_B,D=D,type='r*edf0'))
        	}
        else if(model_A@Model %in% c('OLS') & model_B@Model=='GWR') {
        	pv<-1-pf(T, (model_B@tSS-model_A@tSS)^2/(model_B$tSS2-model_A@tSS2), model_B$edf)
        	        	Accept_HO=(pv<Sidak_cor_MGWR(alpha=0.05,model_A,model_B,D=D,type='r*edf0'))

        	}
        else if(model_A$Model=='MGWRSAR_0_kc_kv' & model_B$Model=='MGWRSAR_0_0_kv') {
        	pv<-1-pf(T, (model_B@tSS-model_A@tSS)^2/(model_B@tSS2-model_A@tSS2), (n-model_B@tSS)^2/(n-model_B@tSS2))
        	        	Accept_HO=(pv<Sidak_cor_MGWR(alpha=0.05,model_A,model_B,D=D,type='r*edf0'))

        	}
        else if(model_A@Model=='OLS' & model_B@Model=='MGWRSAR_0_0_kv') {
        	pv<-1-pf(T, model_A@edf-model_B@edf, model_B@edf)
        	        Accept_HO=(pv<Sidak_cor_MGWR(alpha=0.05,model_A,model_B,D=D,type='r*edf0'))
        	}
        else if(model_A@Model=='SAR' & model_B@Model=='MGWRSAR_0_0_kv') {
        	pv<-1-pf(T, model_A@edf-model_B@edf, model_B@edf)
        	        	Accept_HO=(pv<Sidak_cor_MGWR(alpha=0.05,model_A,model_B,D=D,type='r*edf0'))
        	}
        else if(model_A@Model=='MGWR' & model_B@Model=='MGWRSAR_0_kc_kv') Accept_HO<-(1-pf(T, 1, model_B$edf)<0.5)
        else if(model_A@Model=='GWR' & model_B@Model=='MGWRSAR_0_0_kv') Accept_HO<-(1-pf(T, 1, model_B$edf)<0.5)
        else if(model_A@Model=='OLS' & model_B@Model=='SAR') Accept_HO<-(1-pf(T, 1, model_B$edf)<0.5)

        z=list(T=T,Accept_HO=Accept_HO)
    } else z=T
    z
}
