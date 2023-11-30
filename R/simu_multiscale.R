#' Estimation of linear and local linear model with spatial autocorrelation model (mgwrsar).
#'
#' The simu_multiscale function is designed for simulating a spatially varying
#' coefficient DGP (Data Generating Process) based on formulations proposed by
#' Fotheringam et al. (2017), Gao et al. (2021), or Geniaux (2024).
#'
#' @usage simu_multiscale(n=1000,myseed=1,type='GG2024',b0_constant=FALSE)
#' @param n   An integer number of observations.
#' @param myseed An integer seed used for the simulation.
#' @param type Type of DGP used 'FT2017', 'Gao2021' or 'GG2024', default 'GG2024'.
#' @param b0_constant A boolean parameter indicating whether the intercept term
#' should be spatially varying (TRUE) or not (FALSE).
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
simu_multiscale<-function(n=1000,myseed=1,type='GG2024',b0_constant=FALSE){
  set.seed(myseed)
  x=runif(n)
  y=runif(n)
  coords=cbind(x,y)
  X1=runif(n,-sqrt(3),sqrt(3))
  X2=rnorm(n)
  eps=rnorm(n,0,0.5^2)
  if(b0_constant) Beta1=3 else Beta1=2*(x+y)
  if(type=='FT2017'){
    Beta2=1+(x+y)/12
    Beta3=1+1/324*(36-(6-x/2)^2)*(36-(6-y/2)^2)
    Y=Beta1+Beta2*X1+Beta3*X2+eps
    mydata=data.frame(Y,X1,X2,Beta1,Beta2,Beta3,eps,x,y)
  } else if(type=='Gao2021'){
    Beta2=4 *sin(sqrt(12*(x-0.5)^2+12*(y-0.5)^2))
    Beta3=64 *(x*y) *(1-x)*(1-y)
    Y=Beta1+Beta2*X1+Beta3*X2+eps
    mydata=data.frame(Y,X1,X2,Beta1,Beta2,Beta3,eps,x,y)
  } else {
    Beta2=64 *(x*y) *(1-x)*(1-y)
    Beta3=4*sin(sqrt((6^2*(x-y)^2)))
    Beta4=4*sin(sqrt((6^3*(-x-y)^2)))
    X3=rnorm(n)
    Y=Beta1+Beta2*X1+Beta3*X2+Beta4*X3+eps
    mydata=data.frame(Y,X1,X2,X3,Beta1,Beta2,Beta3,Beta4,eps,x,y)
  }
  mydata$Intercept=1
  list(mydata=mydata,coords=coords)
}
