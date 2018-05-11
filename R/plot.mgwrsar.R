#'  plot_mgwrsar plots the value of local paramaters of a mgwrsar models using a leaflet map.
#' @usage plot_mgwrsar(model,type='coef',var=NULL,SP=NULL,SP_id=NULL,
#' proj=NULL,mypalette= "RdYlGn",opacity=1,fopacity=1,radius=1500)
#' @param model   a mgwsar model.
#' @param type   default 'coef', for plotting the value of the coefficients. Local t-Student could also be plot using 't_coef'.
#' @param var   Names of variable to plot.
#' @param SP   A spdf object.
#' @param SP_id   Id regions for spdf object.
#' @param proj   A CRS projection.
#' @param mypalette   A leaflet palette.
#' @param opacity    Opacity of border color.
#' @param fopacity   Opacity of fill color.
#' @param radius   radius of circle for plot of points.
#' @return A Interactive Web Maps with local parameters plot and Open Street Map layer.
#' @seealso  MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, predict_mgwrsar, kernelW_C
#' @examples
#' \donttest{
#' library(mgwrsar)
#' data(data_mgwrsar)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' W=KNN(coord,4)
#' model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata, coord=coord,
#' fixed_vars=NULL,kernels=c('gauss'),H=0.13,
#' Model = 'GWR', control=list(SE=TRUE))
#' summary_mgwrsar(model_GWR)
#' plot_mgwrsar(model_GWR,type='t_coef',var='X1')
#' }
plot_mgwrsar<-
  function(model,type='coef',var=NULL,SP=NULL,SP_id=NULL,proj=NULL,mypalette= "RdYlGn",opacity=1,fopacity=1,radius=1500) {
if(class(model)!='mgwrsar') stop('A mgwrsar class object is needed')
# x=x_MGWRSAR_0_0_kv
if(type=='t_coef') {
toplot<-data.frame(model$Betav/model$sev)
names(toplot)=colnames(model$Betav)
mytitle=paste("t_Beta(ui,vi)",var)
} else if(type=='residuals') {
toplot<-data.frame(model$residuals)
names(toplot)<-mytitle<-var<-'residuals'
} else if(type=='fitted') {
toplot<-data.frame(model$fit)
names(toplot)<-mytitle<-var<-'fitted'
} else { toplot<-data.frame(model$Betav)
names(toplot)=colnames(model$Betav)
mytitle=paste("Beta(ui,vi)",var)
}
if(!is.null(model$id)) toplot$id=model$id else toplot$id=1:nrow(toplot)
#toplot$id_insee=commune=commune
if(is.null(SP)) {
coordinates(toplot)<-model$coord
if(is.null(proj)) proj4string(toplot)<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" else proj4string(toplot)<-proj
toplot@data$X1=toplot@data[,var]
} else {
SP@data$X1=toplot[match(SP@data[,SP_id],toplot$id),var]
toplot=SP
toplot$id=SP@data[,SP_id]
}
if(proj4string(toplot)!="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") toplot <- spTransform(toplot,CRSobj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# graph carto_pos
palnum <- colorNumeric(palette = mypalette, domain = NULL,reverse=FALSE)
labels <- sprintf(
  paste("<strong>Identifiant</strong>: %s </br><strong>",var,"</strong>: %g"),toplot@data$id,toplot@data$X1
) %>% lapply(htmltools::HTML)

carto <- leaflet(toplot) %>% addTiles()  %>% setView(lng =mean(coordinates(toplot)[,1]), lat = mean(coordinates(toplot)[,2]), zoom = 8)

if(class(toplot)=='SpatialPolygonsDataFrame'){
carto<-carto %>% addPolygons(color = "", weight = 3,opacity = opacity, fillOpacity = fopacity,fillColor = ~palnum(X1),highlightOptions = highlightOptions(color = "white", weight = 3,bringToFront = TRUE),label = labels,labelOptions = labelOptions(style = list("font-weight" = "normal", padding = "3px 8px"),textsize = "15px",direction = "auto"))
} else if(class(toplot)=='SpatialPointsDataFrame'){
carto<-carto %>% addCircles(radius =radius,color = "", weight = 3,opacity = opacity, fillOpacity = fopacity,fillColor = ~palnum(X1),highlightOptions = highlightOptions(color = "white", weight = 3,bringToFront = TRUE),label = labels,labelOptions = labelOptions(style = list("font-weight" = "normal", padding = "3px 8px"),textsize = "15px",direction = "auto"))
}
carto<-carto %>%addLegend(pal = palnum, values = ~X1,bins =11,labFormat=labelFormat(digits=5), opacity = 1, title = mytitle)
carto
}
