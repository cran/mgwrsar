#' prep_w
#' to be documented
#' @usage prep_w(H,kernels,Type='GD',adaptive=FALSE,dists=NULL,indexG=NULL,
#' alpha=1,K=4,theta=0.5)
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param Type  Type of Kernel ('T','GDT','GD'), default 'GD'
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel
#' @param dists  List of precomputed Matrix of spatial distances, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL
#' @param alpha A ratio between 0 and 1 for GDT kernels,default 1
#' @param theta A ratio between 0 and 1 for seasonnal GDT kernels, default 0.5

#' @noRd
#' @return to be documented
prep_w<-function(H,kernels,Type='GD',adaptive=FALSE,dists=NULL,indexG=NULL,alpha=1,K=4,theta=0.5){
  normW<-function(x){
    rs=rowSums(x)
    x[rs>0]<-x[rs>0]/rs[rs>0]
    x
  }
  mycycle=function(x) ifelse(x<=cycling,x,pmin(abs(x-trunc(x/cycling)*cycling),abs(cycling-x%%cycling)))  ### adaptive : if adaptive round(H) and rename kernel
  for(j in 1:length(adaptive)){
    case = adaptive[j]
    if(case) {
      if(Type=='GD' & any(H!=H[1])) H=round(H) else H[j]=round(H[j]) # old condition for local H, to clean
      kernels[j]= paste0(kernels[j],'_adapt_sorted')
    }
  }
  if(Type=='T') {
    if(kernels=='past') Wd<-(dists[['dist_t']]>0)*1 ### only past observations ~ rectangle kernel
    else  {
      Wd=do.call(kernels[1],args=list(dists[['dist_t']],H[1]))
    }
  } else {
    Wd=do.call(kernels[1],args=list(dists[['dist_s']],H[1]))
  }
  #Wd=normW(Wd)
  if(Type=='GDT') {
    kernels_t<-unlist(str_split(kernels[2], '_'))[1]
    format_t<-unlist(str_split(kernels[2], '_'))[2] ## a modifier passer dans control
    cycling<-as.numeric(unlist(str_split(kernels[2], '_'))[3])  ## a modifier passer dans control
    if(adaptive[2]) stop('Only non adptive kernel are allowed for time when type=GDT')
    ## on transforme les distances pour l'adapter au cycle
    if(!is.na(cycling)) {
      same_cycle<-(dists[['dist_t']]<=cycling)*1
      dist_t=mycycle(dists[['dist_t']])
      } else dist_t=dists[['dist_t']]

    wt=do.call(kernels_t,args=list(dist_t,H[2]))
    if(!is.na(cycling)) {wt=(theta*same_cycle+(1-theta)*(1-same_cycle))*wt}
    #wt=normW(wt)
    #### noyau temporel gaussien non adaptatif
    #### spatial + temporel + spatio-temp ?
    #### autre année ?


    ## A ce stade la ponderation est symétrique dans le temps
    if(!is.na(format_t)){
      if(format_t=='past') { ### only past observations for i> H[2]
        past=(dists[['dist_t']]>=0)*1
        id<-which(rowSums(past)>K)
        if(length(id)>0) wt[id,]<-wt[id,]*past[id,] else wt<-wt*past
      }
    }
    #browser()
    wdt<-Wd*wt ## spatio-temporal
   # wdt<-normW(wdt)
    ## si alpha<1 on rajouter un poids egal à tous les observations
    ## A CORRIGER POUR NN != n
    if(alpha<1) {
      wdt<-alpha*wdt + (1-alpha)*1/nrow(wdt)
    }

    if(alpha<1 & !is.na(format_t)){
      if(format_t=='past') {
      past<-normW(past)
      Wd=alpha*wdt+(1-alpha)*past ## spatio-temp + past aspatial
    } else  Wd=wdt
    } else  Wd=wdt
  }
  Wd=Wd/rowSums(Wd)
  list(indexG=indexG,Wd=Wd,dists=dists)
}
