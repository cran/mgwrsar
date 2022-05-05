#' A sampling method for bootstrap test with spatial data using distance bisquare kernel.
#'
#' @usage sample_spat(D,h,focal='median')
#' @param D  A square matrix of distance.
#' @param h  A distance bandwidth.
#' @param focal  If focal is NULL, the ith residuals can not be selected for the ith observation, if focal=='median' it could be selected with a probability equals to the median probability of neighbors.
#' @noRd
#' @return a vector of bootstrapped index
sample_spat <-
function(D,h,focal='median'){
mysample=c()
	for(i in 1:nrow(D)){
	h0=h
	while(sum(bisq(D[i,], h0)>0)<10) h0=h0*1.2
	prob_i= bisq(D[i,], h0)
	#cat('sum(prob_i>0)',sum(prob_i>0))
	if(focal=='median') prob_i[i]<-median(prob_i[-i][prob_i[-i]>0]) else prob_i[i]<-0
	mysample<-c(mysample,sample(1:nrow(coord),1,prob=prob_i))
	}
mysample
}
