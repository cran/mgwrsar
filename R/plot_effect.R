#' plot_effect
#' plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables.
#' @usage plot_effect(model,sampling=TRUE,nsample=2000,title='')
#' @param model a model of mgwrsar class with some spatially varying coefficients.
#' @param sampling Bolean, if nrow(model$Betav)> 5000 a sample of size nsample is randomly selected, default TRUE.
#' @param nsample integer,  size of the sample if sampling is TRUE, default 2000.
#' @param title a title for the plot.
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coord=as.matrix(mydata[,c("x_lat","y_lon")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 8 nearest neighbors with 0 in diagonal
#'  model_GWR0<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#'  fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model = 'GWR',control=list(SE=TRUE))
#'  plot_effect(model_GWR0)
#' }
plot_effect<-function(model,sampling=TRUE,nsample=2000,title=''){
  aes_bx<-list(
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title.x=ggplot2::element_blank(),
      # axis.text.x=element_blank(),axis.ticks.x=element_blank(),
      axis.title.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      # axis.text.y=element_blank(),
      axis.text=ggplot2::element_text(size=10),
      # legend.position = c(0.93, 0.3),
      # legend.direction="vertical",
      legend.position="none"
      # legend.text=element_text(size=8))
    )
  )

if(is.null(model$XV)) stop('No spatially varying coefficients in this model.')
if(sampling & nrow(model$XV)>5000) {
  idx=sample(1:nrow(model$Betav),nsample)
  model$Betav=model$Betav[idx,]
  model$XV=model$XV[idx,]
  model$sev=model$sev[idx,]
}

n=nrow(model$XV)
col_lambda<-which(colnames(model$Betav)=='lambda')
if(length(col_lambda)>0) {
  indexcol=(1:ncol(model$Betav))[-which(colnames(model$Betav)=='lambda')]
  Wy=as.numeric(model$Betav[,col_lambda]* model$W %*%model$Y)
  } else indexcol=(1:ncol(model$Betav))

if(!is.null(model$sev)) {
for(i in 1:ncol(model$Betav)){
  xv=colnames(model$Betav)[i]
  if(i==1) res=data.frame(variable=rep(xv,n),eff=model$XV[,xv]*model$Betav[,xv],signif_95=abs(model$Betav[,xv]/model$sev[,xv])>1.96) else res=rbind(res,data.frame(variable=rep(xv,n),eff=model$XV[,xv]*model$Betav[,xv],signif_95=abs(model$Betav[,xv]/model$sev[,xv])>1.96))
}
if(length(col_lambda)>0) res=rbind(res,data.frame(variable=rep('Wy',n),eff=Wy),signif_95=rep(NA,n))

gAS<-ggplot2::ggplot(res,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff,colour = signif_95))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::scale_colour_manual(values = c("red", 'black'))+ggplot2::geom_vline(xintercept=0)+aes_bx
} else {
  for(i in 1:ncol(model$Betav[,indexcol])){
    xv=colnames(model$Betav)[i]
    if(i==1) res=data.frame(variable=rep(xv,n),eff=model$XV[,xv]*model$Betav[,xv]) else res=rbind(res,data.frame(variable=rep(xv,n),eff=model$XV[,xv]*model$Betav[,xv]))
  }
  if(length(col_lambda)>0) res=rbind(res,data.frame(variable=rep('Wy',n),eff=Wy))
  gAS<-ggplot2::ggplot(res,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::geom_vline(xintercept=0)+aes_bx
}
gAS+ ggplot2::ggtitle(title)
}
