#' built_Vseq
#' to be documented
#' @usage built_Vseq(d,h)
#' @param env an environment
#' @return an environment.
#' @noRd
#'
built_Vseq<-function(env=parent.frame()){
with(env,{
  if(is.null(first_nn)) first_nn=n
  first_nn<-min(first_nn,control$NN)
  if(is.null(V)){
    if(type=='proportional'){
      alpha2=exp(log(minv/first_nn)/(nns+1))
      V=sapply(1:(nns+1),function(x) round(first_nn*(alpha2)^x))
    }  else {
      v=round(first_nn/nns)
      V=seq(first_nn,minv,by=-v)
    }
    V<-unique(V[V>=minv])
  } else {
    m=V[1]
    V<-V[V<(n-2)]
  }
  if(!is.null(minv))  V=V[V>minv]
  if(control$adaptive[1]) max_dist=n else {
    max_dist=max(G$dists[['dist_s']])
    if(is.null(min_dist)) min_dist= quantile(G$dists[['dist_s']][,2],0.9)
  }
  if(control$Type=='GDT'){
    cycling<-as.numeric(unlist(str_split(kernels[2], '_'))[3])
    if(length(cycling)>0){
    Vt=quantile(abs(G$dists[['dist_t']])[abs(G$dists[['dist_t']])<cycling],seq(1,0.05,by=-0.05))
    Vt<-c(quantile(abs(G$dists[['dist_t']]),0.95),Vt)
    } else Vt=quantile(abs(G$dists[['dist_t']]),seq(1,0.05,by=-0.05))
    #max_dist=quantile(G$dists[['dist_s']],0.8)
    #alpha=exp(log(min_dist/max_dist)/(nns+1))
    #V=sapply(1:(nns+1),function(x) round(max_dist*(alpha)^x,digits=10))
  }
  if(!control$adaptive[1]){
    #browser()
    V=sapply(V,function(x) median(G$dists[['dist_s']][,x]))
    #max_dist=quantile(G$dists[['dist_s']],0.8)
    #alpha=exp(log(min_dist/max_dist)/(nns+1))
    #V=sapply(1:(nns+1),function(x) round(max_dist*(alpha)^x,digits=10))
  }


  min_dist=min(V)
  l=length(V)
  if(Model!='atds_gwr') V=c(max_dist,V)
 })
}





