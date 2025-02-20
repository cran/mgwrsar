#' plot_effect
#' plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables.
#' @usage plot_effect(model,sampling=TRUE,nsample=2000,nsample_max=5000,title='')
#' @param model a model of mgwrsar class with some spatially varying coefficients.
#' @param sampling Bolean, if nrow(model@Betav)> nsample_max a sample of size nsample is randomly selected, default TRUE.
#' @param nsample integer,  size of the sample if sampling is TRUE, default 2000.
#' @param nsample_max integer, size max to engage sampling if sampling is TRUE, default 5000.
#' @param title a title for the plot.
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 8 nearest neighbors with 0 in diagonal
#'  model_GWR0<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coords=coords,
#'  fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model = 'GWR',control=list(SE=TRUE))
#'  plot_effect(model_GWR0)
#' }
plot_effect<-function(model,sampling=TRUE,nsample=2000,nsample_max=5000,title=''){
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

  if(is.null(model@XV)) stop('No spatially varying coefficients in this model.')
  n=nrow(model@XV)
  m=ncol(model@XV)
  #cat('length(model@Betac)',length(model@Betac))
  if(length(model@Betac)>0) {
    model@XV= cbind(model@XV,model@XC)
    model@Betav=cbind(model@Betav,matrix(model@Betac,nrow=n,ncol=length(model@Betac),byrow =TRUE))
    colnames(model@Betav)[-(1:m)]<-names(model@Betac)
  }
  if(!is.null(model@sev)) {
    if(length(model@Betac)>0) {
      model@sev=cbind(model@sev,matrix(model@se,nrow=n,ncol=length(model@se),byrow =TRUE))
      colnames(model@sev)[-(1:m)]<-names(model@Betac)
    }
  }

  if(sampling & nrow(model@XV)>nsample_max) {
    idx=sample(1:nrow(model@Betav),nsample)
    model@Betav=model@Betav[idx,]
    model@XV=model@XV[idx,]
    model@sev=model@sev[idx,]
  }
  n=nrow(model@XV)

  col_lambda<-which(colnames(model@Betav)=='lambda')
  if(length(col_lambda)>0) {
    indexcol=(1:ncol(model@Betav))[-which(colnames(model@Betav)=='lambda')]
    Wy=as.numeric(model@Betav[,col_lambda]* model@W %*%model@Y)
  } else indexcol=(1:ncol(model@Betav))

  if(!is.null(model@sev)) {
    for(i in 1:ncol(model@Betav)){
      xv=colnames(model@Betav)[i]
      if(i==1) res=data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],signif_95=abs(model@Betav[,xv]/model@sev[,xv])>1.96,spv=TRUE) else if(i<=m) res=rbind(res,data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],signif_95=abs(model@Betav[,xv]/model@sev[,xv])>1.96,spv=TRUE)) else res=rbind(res,data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],signif_95=abs(model@Betav[,xv]/model@sev[,xv])>1.96,spv=FALSE))
    }
    if(length(col_lambda)>0) res=rbind(res,data.frame(variable=rep('Wy',n),eff=Wy),signif_95=rep(NA,n))

    res$ord<-reorder(res$variable, res$eff, FUN = function(x) abs(max(x)-min(x)))

    res_v<-res %>% dplyr::filter(spv==TRUE) %>% dplyr::arrange(desc(ord),desc(signif_95))
    if(res_v$signif_95[1]) col_v=c( 'black',"red") else col_v=c("red",'black')
    if(length(model@Betac)>0) {
      res_c<-res %>% dplyr::filter(spv==FALSE) %>% dplyr::arrange(desc(ord),desc(signif_95))
      if(res_c$signif_95[1]!=res_v$signif_95[1]) col_c=col_v else col_c=col_v[c(2,1)]
    }


    gAS_v<-ggplot2::ggplot(res_v,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff,colour = signif_95))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::scale_colour_manual(values = col_v)+ggplot2::geom_vline(xintercept=0)+aes_bx+ ggplot2::ggtitle('Varying coefficients')



    if(length(model@Betac)>0) gAS_c<-ggplot2::ggplot(res_c,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff,colour = signif_95))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::scale_colour_manual(values = col_c)+ggplot2::geom_vline(xintercept=0)+aes_bx + ggplot2::ggtitle('Constant coefficients')


  } else {
    for(i in 1:ncol(model@Betav[,indexcol])){
      xv=colnames(model@Betav)[i]
      if(i==1) res=data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],spv=TRUE) else if(i<=m) res=rbind(res,data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],spv=TRUE)) else res=rbind(res,data.frame(variable=rep(xv,n),eff=model@XV[,xv]*model@Betav[,xv],spv=FALSE))
    }
    if(length(col_lambda)>0) res=rbind(res,data.frame(variable=rep('Wy',n),eff=Wy))

    res_v<-res %>% filter(spv==TRUE) %>% arrange(desc(ord),desc(signif_95))
    if(length(model@Betac)>0) res_c<-res %>% filter(spv==FALSE) %>% arrange(desc(ord),desc(signif_95))

    gAS_v<-ggplot2::ggplot(res_v,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::geom_vline(xintercept=0)+aes_bx + ggplot2::ggtitle('Varying coefficients')

    if(length(model@Betac)>0) gAS_c<-ggplot2::ggplot(res_c,ggplot2::aes(y=reorder(variable, eff, FUN = function(x) abs(max(x)-min(x))),x=eff))+ggplot2::geom_jitter(width = 0.01,size=0.2)+ggplot2::geom_vline(xintercept=0)+aes_bx + ggplot2::ggtitle('Constant coefficients')


  }
  if(length(model@Betac)>0){
    xrange<-ggplot_build(gAS_v)$layout$panel_scales_x[[1]]$range$range
    gAS_c<-gAS_c+xlim(xrange[1],xrange[2])
    grid.arrange(gAS_v,gAS_c,ncol=2,top = textGrob(title,gp=gpar(fontsize=20,font=3)))
  } else  gAS_v
}
