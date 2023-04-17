#'  plot_mgwrsar plots the value of local paramaters of a mgwrsar models using a leaflet map.
#' @usage plot_mgwrsar(model,type='coef',var=NULL,crs=NULL,mypalette= "RdYlGn",opacity=0.5
#' ,fopacity=0.5,nbins=8,radius=500,mytile='Stamen.TonerBackground',myzoom=8,
#' myresolution=150,LayersControl=TRUE,myzoomControl=TRUE,mytile2=NULL,ScaleBar=NULL,
#' ScaleBarOptions=list(maxWidth = 200, metric = TRUE,imperial = FALSE,
#' updateWhenIdle = TRUE),MyLegendTitle=NULL,lopacity=0.5)
#' @param model   a mgwsar model.
#' @param type   default 'coef', for plotting the value of the coefficients. Local t-Student could also be plot using 't_coef', residuals using 'residuals' and fitted using 'fitted'.
#' @param var   Names of variable to plot.
#' @param crs   A CRS projection.
#' @param mypalette   A leaflet palette.
#' @param opacity    Opacity of border color.
#' @param fopacity   Opacity of fill color.
#' @param radius   radius of circle for plot of points.
#' @param nbins nbins.
#' @param mytile tile 1.
#' @param myzoom level of zoom for tile 1.
#' @param myresolution resolution for tile 1.
#' @param LayersControl layers contols.
#' @param myzoomControl zoem control.
#' @param mytile2  tile 2.
#' @param ScaleBar ScaleBar.
#' @param ScaleBarOptions options for ScaleBar.
#' @param MyLegendTitle Legend title.
#' @param lopacity opacity for legend.
#' @return A Interactive Web Maps with local parameters plot and Open Street Map layer.
#' @seealso  MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, predict_mgwrsar, kernel_matW
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coord=as.matrix(mydata[,c("x_lat","y_lon")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 4 nearest neighbors with 0 in diagonal
#'  model_GWR0<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#'  fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model='GWR',control=list(SE=TRUE))
#'  summary_mgwrsar(model_GWR0)
#'  plot_mgwrsar(model_GWR0,type='B_coef',var='X2')
#'  plot_mgwrsar(model_GWR0,type='t_coef',var='X2')
#' }
plot_mgwrsar<-
  function(model,type='coef',var=NULL,crs=NULL,mypalette= "RdYlGn",opacity=0.5,fopacity=0.5,nbins=8,radius=500,mytile='Stamen.TonerBackground',myzoom=8,myresolution=150,LayersControl=TRUE,myzoomControl=TRUE,mytile2=NULL,ScaleBar=NULL,ScaleBarOptions=list(maxWidth = 200, metric = TRUE,imperial = FALSE, updateWhenIdle = TRUE),MyLegendTitle=NULL,lopacity=0.5) {
if(!is(model,'mgwrsar')) stop('A mgwrsar class object is needed')
# x=x_MGWRSAR_0_0_kv
if(type=='t_coef') {
toplot<-data.frame(model$Betav[,var]/model$sev[,var])
names(toplot)='X1'
mytitle=paste("t_Beta(ui,vi)",var)
} else if(type=='residuals') {
toplot<-data.frame(model$residuals)
names(toplot)<-'X1'
mytitle<-var<-'residuals'
} else if(type=='fitted') {
toplot<-data.frame(model$fit)
names(toplot)<-'X1'
mytitle<-var<-'fitted'
} else { toplot<-data.frame(model$Betav[,var])
names(toplot)='X1'
mytitle=paste("Beta(ui,vi)",var)
}
if(!is.null(MyLegendTitle)) mytitle=MyLegendTitle
if(!is.null(model$id)) toplot$id=model$id else toplot$id=1:nrow(toplot)
#toplot$id_insee=commune=commune
toplot = st_as_sf(cbind(toplot,model$coord), coords = c('x_lat','y_lon'), remove = TRUE)
if(is.null(crs)) toplot <- st_set_crs(toplot, '+proj=longlat +datum=WGS84')  else toplot <- st_transform(st_set_crs(toplot,crs),'+proj=longlat +datum=WGS84')
# graph carto_pos
palnum <- colorQuantile(palette = mypalette, domain =NULL,n=nbins,reverse = TRUE)
labels <- paste0("<strong>Identifiant</strong>: ",toplot$id," </br><strong>",var,"</strong>: " ,as.numeric(unlist(data.frame(toplot$X1)))) %>% lapply(htmltools::HTML)
carto <- leaflet(toplot,options = leafletOptions(resolutions = myresolution,zoomSnap = 0.25, zoomDelta=0.25,zoomControl = myzoomControl)) %>% addMapPane(name = "points", zIndex = 410)  %>% addMapPane(name = "label", zIndex = 420) %>% setView(lng =(min(st_coordinates(toplot)[,1])+max(st_coordinates(toplot)[,1]))/2, lat = (min(st_coordinates(toplot)[,2])+max(st_coordinates(toplot)[,2]))/2, zoom = myzoom) %>% addProviderTiles(mytile,group = mytile)
carto<-carto %>% addCircles(radius =radius,color = "", weight = 3,opacity = opacity, fillOpacity = fopacity,fillColor = ~palnum(X1),highlightOptions = highlightOptions(color = "white", weight = 3,bringToFront = TRUE),label = labels,labelOptions = labelOptions(style = list("font-weight" = "normal", padding = "3px 8px"),textsize = "15px",direction = "auto"),group = "my points",options = leafletOptions(pane = "points",resolutions = myresolution,zoomSnap = 0.25, zoomDelta=0.25,zoomControl = myzoomControl))
if(!is.null(mytile2)) carto <- carto %>% addProviderTiles(mytile2,group = "map labels",options = leafletOptions(pane = "label"))
if(LayersControl) carto <- carto %>%addLayersControl(baseGroups = mytile, overlayGroups = c("map labels", "my points"))

if(!is.null(ScaleBar)) carto <- carto %>%addScaleBar(position = ScaleBar,options = ScaleBarOptions)

carto<-carto %>%addLegend(pal = palnum, values = ~X1,labFormat = function(type, cuts, p) {
  n = length(cuts)
  p = paste0(round(p * 100), '%')
  cuts = paste0(formatC(cuts[-n]), " - ", formatC(cuts[-1]))
  # mouse over the legend labels to see the percentile ranges
  paste0(
    '<span title="', p[-n], " - ", p[-1], '">', cuts,
    '</span>'
  )
}, opacity = lopacity, title = mytitle)
carto
}
