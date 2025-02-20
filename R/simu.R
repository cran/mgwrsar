#' Estimation of linear and local linear model with spatial
#' autocorrelation model (mgwrsar).
#'
#' The simu_multiscale function is designed for simulating a spatially varying
#' coefficient DGP (Data Generating Process) based on formulations proposed by
#' Fotheringam et al. (2017), Gao et al. (2021), or Geniaux (2024).
#'
#' @usage simu_multiscale(n=1000,myseed=1,type='GG2024',constant=NULL,
#' nuls=NULL,config_beta='default',config_snr=0.7,config_eps='normal',
#' ratiotime=1)
#' @param n   An integer number of observations
#' @param myseed An integer seed used for the simulation.
#' @param type Type of DGP used 'FT2017', 'Gao2021' or 'GG2024', default 'GG2024'.
#' @param constant A boolean parameter indicating whether the intercept term
#' should be spatially varying (TRUE) or not (FALSE).
#' @param nuls A vector of null parameters, default NULL
#' @param config_beta name of the type of spatial pattern of Beta coefficients
#' @param config_snr a value of signal noise ratio
#' @param config_eps name of the distribution of error ('normal','unif' or 'Chi2')
#' @param ratiotime multiplicating factor, for spacetime DGP.
#' @return A named list with simulated data ('mydata') and coords ('coords')
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  library(ggplot2)
#'  library(gridExtra)
#'  library(grid)
#'  simu=simu_multiscale(1000)
#'  mydata=simu$mydata
#'  coords=simu$coords
#'  p1<-ggplot(mydata,aes(x,y,col=Beta1))+geom_point() +scale_color_viridis_c()
#'  p2<-ggplot(mydata,aes(x,y,col=Beta2))+geom_point() +scale_color_viridis_c()
#'  p3<-ggplot(mydata,aes(x,y,col=Beta3))+geom_point() +scale_color_viridis_c()
#'  p4<-ggplot(mydata,aes(x,y,col=Beta4))+geom_point() +scale_color_viridis_c()
#'  grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2, top = textGrob("DGP Geniaux (2024)"
#'  ,gp=gpar(fontsize=20,font=3)))
#' }
simu_multiscale<-function(n=1000,myseed=1,type='GG2024',constant=NULL,nuls=NULL,config_beta='default',config_snr=0.7,config_eps='normal',ratiotime=1){
  set.seed(myseed)
  get_snr<-function(k,config_snr,XB,config_eps){
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif') eps=k*runif(n, min = -3, max = 3)
    else if(config_eps=='Chi2') eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    abs(config_snr-(sum(XB^2))/sum((XB+eps)^2))
  }

  if(config_beta=='All_B3') type='All_B3'
  if(config_beta=='B4_nonconstant') constant=c(1,2,3)
  if(config_beta=='B1_constant') constant=1
  if(config_beta=='B2_constant') constant=2
  if(config_beta=='B3_constant') constant=3
  if(config_beta=='B4_constant') constant=4
  if(config_beta=='B1_nul') nuls=1
  if(config_beta=='B2_nul') nuls=2
  if(config_beta=='B3_nul') nuls=3
  if(config_beta=='B4_nul') nuls=4
  if(config_beta=='GG2024_univ') type='GG2024_univ'
  if(config_beta=='low_h') type='GG2024_low'
  if(config_beta=='high_h') type='GG2024_high'
  if(config_beta=='FT2017') type='FT2017'
  if(config_beta=='Gao2021') type='Gao2021'
  if(config_beta=='GG2024_2') type='GG2024_2'
  if(config_beta=='GG2024_8') type='GG2024_8'
  if(config_beta=='GG2024_16') type='GG2024_16'
  if(config_beta=='GG2024_32') type='GG2024_32'
  if(config_beta=='spatiotemp') type='spatiotemp'


  x=runif(n)
  y=runif(n)
  coords=cbind(x,y)
  X1=runif(n,-sqrt(3),sqrt(3))
  X2=rnorm(n)

  Beta1=3*(x+y)-1
  if(type=='FT2017'){
    Beta1=3
    Beta2=1+(x*25+y*25)/12
    Beta3=1+1/324*(36-(6-x*25/2)^2)*(36-(6-y*25/2)^2)

  } else if(type=='Gao2021'){
    Beta2=4 *sin(sqrt(12*(x-0.5)^2+12*(y-0.5)^2))
    Beta3=84 *(x*y) *(1-x)*(1-y)

  }  else if(type=='GG2024_univ' ) {
    Beta1=3*x-1
    Beta2=84 *(x/2) *(1-x)*0.5
    Beta3=4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    Beta4=4*sin(sqrt((6^3*(-x-0.25)^2)))+2*(x-0.5)
    X3=rnorm(n)

  } else if(type=='GG2024_low' ) {
    Beta1=3*x-1
    Beta2=84 *(x*y) *(1-x)*(1-y)
    Beta3=3*(x-y)-1#4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    Beta4=3*(y-x)-1
    X3=rnorm(n)

  } else if(type=='GG2024_high' ) {
    Beta1=4*x*sin(sqrt((6^2*(-0.4+x+y)^4)))
    Beta2=84 *(x*y) *(1-x)*(1-y)
    Beta3=4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    Beta4=4*sin(sqrt((6^3*(-x-y)^2)))+2*(x-0.5)
    X3=rnorm(n)

  }  else if(type=='All_B3' ) {
    Beta1<-3
    Beta2<-Beta4<-Beta3<-4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    X3=rnorm(n)
  }  else if(type=='spatiotemp' ) {
    time=seq(0,1, length.out = n)
    timeratio=1+time*ratiotime
    Beta1=(3*(x+y)-1)*timeratio
    Beta2=(84 *(x*y) *(1-x)*(1-y))*timeratio
    Beta3=4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    Beta4=4*sin(sqrt((6^3*(-x-y)^2)))+2*(x-0.5)
    X3=rnorm(n)
  } else {
    Beta2=84 *(x*y) *(1-x)*(1-y)
    Beta3=4*x*sin(sqrt((6^2*(-0.4+x*2)^4)))
    Beta4=4*sin(sqrt((6^3*(-x-y)^2)))+2*(x-0.5)
    X3=rnorm(n)
  }
  if(length(constant)==1){
    if(constant==1) Beta1=3
    if(constant==2) {
      Beta2=3;
    }
    if(constant==3)  {
      Beta3=3;
    }
    if(constant==4)  {
      Beta4=3;
    }
  } else {
    for(i in 1:length(constant)) assign(paste0('Beta',constant[i]),3)
  }

  if(length(nuls)==1 & length(constant)==0){
    if(nuls==1) Beta1=0
    if(nuls==2) {
      Beta2=0;
    }
    if(nuls==3)  {
      Beta3=0;
    }
    if(nuls==4)  {
      Beta4=0;
    }
  } else {
    for(i in 1:length(nuls)) assign(paste0('Beta',nuls[i]),0)
  }

  if(type %in% c('FT2017','Gao2021')){
    XB=Beta1+Beta2*X1+Beta3*X2
    k=(optimize(get_snr,lower=1,upper=100,config_snr=config_snr,XB=XB,config_eps=config_eps))$minimum
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif')  eps=k*runif(n, min = -3, max = 3) else if(config_eps=='Chi2')  eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    Y=Beta1+Beta2*X1+Beta3*X2+eps
    mydata=data.frame(Y,X1,X2,Beta1,Beta2,Beta3,eps,x,y)
  } else if(type %in% c('GG2024_univ')){
    if(config_snr>=0.9) eps=rnorm(n,0,0.1) else if(config_snr>0.7) eps=rnorm(n,0,0.25) else if(config_snr>0.5) eps=rnorm(n,0,0.5)
    Y=Beta1+Beta2*X1+Beta3*X2+Beta4*X3+eps
    Y1=Beta1+eps
    Y2=Beta2*X1+eps
    Y3=Beta3*X2+eps
    Y4=Beta4*X3+eps

    mydata=data.frame(Y,Y1,Y2,Y3,Y4,X1,X2,X3,Beta1,Beta2,Beta3,Beta4,eps,x,y)
  } else if(type=='GG2024_2'){
    XB=Beta1+Beta2*X1
    k=(optimize(get_snr,lower=1,upper=100,config_snr=config_snr,XB=XB,config_eps=config_eps))$minimum
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif')  eps=k*runif(n, min = -3, max = 3) else if(config_eps=='Chi2')  eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    Y=XB+eps
    mydata=data.frame(Y,X1,Beta1,Beta2,eps,x,y)
  }  else if(type=='GG2024_8'){
    X4=runif(n,-sqrt(3),sqrt(3))
    X5=rnorm(n)
    X6=runif(n,-sqrt(3),sqrt(3))
    X7=rnorm(n)
    Beta5=1+(x*25+y*25)/12
    Beta6=1+1/324*(36-(6-x*25/2)^2)*(36-(6-y*25/2)^2)
    Beta7=4 *sin(sqrt(12*(x-0.5)^2+12*(y-0.5)^2))
    Beta8=84 *(x*y) *(1-x)*(1-y)
    XB=Beta1+Beta2*X1+Beta3*X2+Beta4*X3+Beta5*X4+Beta6*X5+Beta7*X6+Beta8*X7
    k=(optimize(get_snr,lower=1,upper=100,config_snr=config_snr,XB=XB,config_eps=config_eps))$minimum
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif')  eps=k*runif(n, min = -3, max = 3) else if(config_eps=='Chi2')  eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    Y=XB+eps
    mydata=data.frame(Y,X1,X2,X3,X4,X5,X6,X7,Beta1,Beta2,Beta3,Beta4,Beta5,Beta6,Beta7,Beta8,eps,x,y)
  }  else if(type=='GG2024_16'){
    X4=runif(n,-sqrt(3),sqrt(3))
    X5=rnorm(n)
    X6=runif(n,-sqrt(3),sqrt(3))
    X7=rnorm(n)
    X8=runif(n,-sqrt(3),sqrt(3))
    X9=rnorm(n)
    X10=runif(n,-sqrt(3),sqrt(3))
    X11=rnorm(n)
    X12=runif(n,-sqrt(3),sqrt(3))
    X13=rnorm(n)
    X14=runif(n,-sqrt(3),sqrt(3))
    X15=rnorm(n)
    Beta5=1+(x*25+y*25)/12
    Beta6=1+1/324*(36-(6-x*25/2)^2)*(36-(6-y*25/2)^2)
    Beta7=4 *sin(sqrt(12*(x-0.5)^2+12*(y-0.5)^2))
    Beta8=84 *(x*y) *(1-x)*(1-y)
    Beta9<-Beta1
    Beta10<-Beta2
    Beta11<-Beta3
    Beta12<-Beta4
    #  A supprimer
    Beta12<-Beta4<-3
    #
    Beta13<-Beta5
    Beta14<-Beta6
    Beta15<-Beta7
    Beta16<-Beta8
    XB=Beta1+Beta2*X1+Beta3*X2+Beta4*X3+Beta5*X4+Beta6*X5+Beta7*X6+Beta8*X7+Beta9*X8+Beta10*X9+Beta11*X10+Beta12*X11+Beta13*X12+Beta14*X13+Beta15*X14+Beta16*X15
    k=(optimize(get_snr,lower=1,upper=100,config_snr=config_snr,XB=XB,config_eps=config_eps))$minimum
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif')  eps=k*runif(n, min = -3, max = 3) else if(config_eps=='Chi2')  eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    Y=XB+eps
    mydata=data.frame(Y,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,Beta1,Beta2,Beta3,Beta4,Beta5,Beta6,Beta7,Beta8,Beta9,Beta10,Beta11,Beta12,Beta13,Beta14,Beta15,Beta16,eps,x,y)
  }  else {
    XB=Beta1+Beta2*X1+Beta3*X2+Beta4*X3
    k=(optimize(get_snr,lower=1,upper=100,config_snr=config_snr,XB=XB,config_eps=config_eps))$minimum
    if(config_eps=='normal') eps=k*rnorm(n,0,1) else if(config_eps=='unif')  eps=k*runif(n, min = -3, max = 3) else if(config_eps=='Chi2')  eps=k*(0.5*rchisq(n, df=2, ncp = 0)-1)
    Y=Beta1+Beta2*X1+Beta3*X2+Beta4*X3+eps
    mydata=data.frame(Y,X1,X2,X3,Beta1,Beta2,Beta3,Beta4,eps,x,y)
  }
  ## ici
  mydata$Intercept=1
  if(config_beta=='spatiotemp') mydata$time=time
  list(mydata=mydata,coords=coords)
}
