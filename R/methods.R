#' Class of mgwrsar Model.
#'
#' @slot Betav matrix, the estimated varying coefficients, dim(n,kv).
#' @slot Betac numeric, the estimated constant coefficients, length kc.
#' @slot Model character, The type of model.
#' @slot fixed_vars  character, a vector with name of constant covarariate.
#' @slot Y  numeric, the dependent variable.
#' @slot XC  matrix, the explanatory variables with constant coefficients.
#' @slot XV  matrix, the explanatory variables with varying coefficients.
#' @slot X  matrix, the explanatory variables.
#' @slot W  SparseMatrix, the spatial weight matrix for spatial dependence.
#' @slot isgcv  logical, if gcv has been computed.
#' @slot edf  numeric, the estimated degrees of freedom.
#' @slot formula  \code{formula}
#' @slot data  dataframe, The dataframe used for computation.
#' @slot Method  character, the estimation technique for computing the models with Spatial Dependence. '2SLS' or 'B2SLS', default '2SLS'.
#' @slot coords matrix, the spatial coordinates of observations.
#' @slot H  numeric, the bandwidth vector.
#' @slot H2  numeric, the time bandwidth vector.
#' @slot kernels  character, the type of kernel.
#' @slot adaptive logical, adaptive kernel.
#' @slot Type character, the type of General Kernel Product.
#' @slot TP numeric, index of target points.
#' @slot SSRtp  numeric, the sum of square residuals for TP.
#' @slot SSR  numeric, the sum of square residuals.
#' @slot residuals  numeric, the vector of residuals.
#' @slot fit  numeric, the vector of fitted values.
#' @slot pred numeric, the vector of predicted values.
#' @slot sev  matrix, local standard error of varying coefficients.
#' @slot se  numeric, standard error of constant coefficients.
#' @slot tS numeric, Trace(S).
#' @slot Shat, hat matrix
#' @slot R_k, list of hat matrix by var
#' @slot h_w numeric, the bandwidth value for computing W, default 0.
#' @slot kernel_w the type of kernel for computing W, default NULL.
#' @slot RMSE numeric, Root Mean Square Error for Target Points.
#' @slot RMSEtp numeric, Root Mean Square Error for all Points.
#' @slot CV numeric, Leave One Out CV.
#' @slot AIC numeric,  Akaike Criteria.
#' @slot AICc numeric, Corrected Akaike Criteria.
#' @slot AICctp numeric, Corrected Akaike Criteria for TP
#' @slot BIC numeric, Bayesian Information Criteria.
#' @slot R2 numeric, R2.
#' @slot R2_adj numeric, adjusted R2.
#' @slot get_ts logical, if trace of hat matrix Tr(S) should be stored.
#' @slot NN  numeric, the maximum number of neighbors for weights computation
#' @slot doMC logical, parallel computation.
#' @slot ncore numeric, number of cores.
#' @slot mycall a call, the call of the model.
#' @slot ctime numeric, the computing times in seconds.
#' @slot HRMSE matrix, RMSE log.
#' @slot HBETA list, estimated BETA at each iteration.
#' @slot loglik numeric, value of loglik.
#' @slot G list, list of neighboring index and distances (knn object from nabor package).
#' @slot V numeric, neighbors sequence for TDS.
#' @slot Vt numeric, neighbors sequence for TDS.
#' @slot Z numeric, time for GDT kernel type
#' @slot TS numeric, Diagonal of Hat Matrix
#' @slot alpha numeric, ratio for GDT kernels
#' @slot theta numeric, ratio for GDT kernels
#
#' @export
setClass("mgwrsar",
         slots= list(
          Betav= "matrix",
          Betac ="numeric",
          Model= "character",
          fixed_vars=  "character",
          Y=  "numeric",
          XC = "matrix",
          XV = "matrix",
          X = "matrix",
          W = "Matrix",
          isgcv=  "logical",
          edf = "numeric",
          formula = "formula",
          data=  "data.frame",
          Method=  "character",
          coords= "matrix",
          H = "numeric",
          H2 = "numeric",
          kernels  ="character",
          h_w="numeric",
          kernel_w  ="character",
          adaptive ="logical",
          Type="character",
          Shat= "matrix",
          R_k= "list",
          TP= "numeric",
          SSRtp=  "numeric",
          SSR=  "numeric",
          residuals=  "numeric",
          fit=  "numeric",
          pred= "numeric",
          sev=  "matrix",
          se=  "numeric",
          tS="numeric",
          TS="numeric",
          RMSE= "numeric",
          RMSEtp= "numeric",
          CV= "numeric",
          AICc= "numeric",
          AICctp= "numeric",
          AIC= "numeric",
          BIC= "numeric",
          R2= "numeric",
          R2_adj= "numeric",
          get_ts= "logical",
          NN=  "numeric",
          doMC ="logical",
          ncore= "numeric",
          mycall= "call",
          ctime= "numeric",
          HRMSE="matrix",
          HBETA="list",
          G="list",
          loglik="numeric",
          V="numeric",
          Vt="numeric",
          Z="numeric",
          alpha="numeric",
          theta="numeric"
          )
)


#' coef for mgwrsar model
#'
#' @param object A model of class \code{\link{mgwrsar-class}}.
#' @param ... coef parameters forwarded.
#' @return A named list with a matrix of varying coefficients and a vector or non varying coefficients.
#' @export
#' @rdname coef.mgwrsar
setMethod("coef",'mgwrsar', function(object,...)
{
list(Betav=object@Betav,Betac=object@Betac)
}
)


#' fitted for mgwrsar model
#'
#' @param object A model of class \code{\link{mgwrsar-class}}.
#' @param ... fitted parameters forwarded.
#' @return A vector of fitted values.
#' @export
#' @rdname fitted.mgwrsar
setMethod("fitted",'mgwrsar', function(object,...)
{
  object@fit
}
)

#' residuals for mgwrsar model
#'
#' @param object A model of class \code{\link{mgwrsar-class}}.
#' @param ... residuals parameters forwarded.
#' @return A vector of residuals.
#' @export
#' @rdname residuals.mgwrsar
setMethod("residuals",'mgwrsar', function(object,...)
{
  object@residuals
}
)


#' summary for mgwrsar model
#'
#' @param object A model of class \code{\link{mgwrsar-class}}.
#' @param ... summary parameters forwarded.
#' @return A summary object.
#' @export
#' @rdname summary.mgwrsar
setMethod("summary",'mgwrsar', function(object,...)
{
  model=object
  if(!is(model,'mgwrsar')) stop("not a mgwrsar object")
  cat("Call:\n")
  print(model@mycall)
  n <- length(model@Y)
  cat("Model:", model@Model, "\n")
  if(!(model@Model %in% c('OLS','GWR','MGWR'))) cat("Method for spatial autocorrelation:", model@Method, "\n")
  cat("Kernels function:", model@kernels, "\n")
  cat("Kernels adaptive:", ifelse(model@adaptive,'YES','NO'), "\n")
  cat("Kernels type:", model@Type, "\n")
  cat("Bandwidth:", model@H, "\n")
  cat("Computation time:", model@ctime, "\n")
  cat("Use of parallel computing:", model@doMC, " ncore =",model@ncore,"\n")
  cat("Use of rough kernel:",ifelse(model@NN<n & (!model@adaptive | (model@adaptive & model@kernels=='gauss')) ,paste0('YES, ',model@NN,' neighbors / ',n),'NO'),"\n")
  cat("Use of Target Points:", ifelse(is.null(model@TP) | length(model@TP)==n,'NO','YES'), "\n")
  if(length(model@TP)!=n) {
    # cat(paste0("Use of ",ifelse(length(model@TP)>1,length(model@K),'single ')," pass of Target Points \n"))
    # cat("Use of iterative mode of TP choice : ",ifelse(model@TPboost,'yes','no')," \n")
    cat("Number of Target Points ", length(unlist(model@TP)), "\n")
  }
  cat("Number of data points:", n, "\n")
  if(length(model@XC)>0) {cat("Constant parameters:", names(model@Betac), "\n")
    cat(as.numeric(model@Betac),'\n')}
  if(length(model@XV)>0){ cat("Varying parameters:", colnames(model@XV), "\n")
    CM <- apply(model@Betav, 2, summary)
    printCoefmat(CM)
  }
  if (length( model@edf)>0) {
    cat("Effective degrees of freedom:", model@edf, "\n")
  }
  if (length(model@tS)>0) cat("AICc:", model@AICc, "\n")

  if (length(model@CV)>0) cat("LOOCV:", model@CV, "\n")

  cat("Residual sum of squares:", model@SSR, "\n")
  cat("RMSE:", model@RMSE, "\n")

  invisible(model)
}
)

#' predict method for mgwrsar model
#'
#' @param object   A model of class \code{\link{mgwrsar-class}}.
#' @param newdata a matrix or data.frame of new data.
#' @param newdata_coords  a matrix of new coordinates, and eventually other variables if a General Kernel Product is used.
#' @param W the spatial weight matrix for models with  spatial autocorrelation.
#' @param type Type for BLUP estimator, default "BPN". If NULL use predictions without spatial bias correction.
#' @param  h_w A bandwidth value for the spatial weight matrix
#' @param kernel_w kernel type for the spatial weight matrix. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"),
#' epanechnikov ("epane"), gaussian ("gauss")) .
#' @param maxobs  maximum number of observations for exact calculation of solve(I- rho*W), default maxobs=4000.
#' @param beta_proj A boolean, if TRUE the function then return a two elements list(Y_predicted,Beta_proj_out)
#' @param method_pred If method_pred = 'TP' (default) prediction is done by recomputing a MGWRSAR model
#' with new-data as target points, else if method_pred in ('tWtp_model','model','shepard') a matrix
#' for projecting estimated betas is used (see details).
#' @param k_extra number of neighboors for local parameter extrapolation if shepard kernel is used, default 8.
#' @param ... predict parameters forwarded.
#' @return A vector of predictions if beta_proj is FALSE or a list with a vector named Y_predicted and a matrix named Beta_proj_out.
#' @details if method_pred ='tWtp_model',  the weighting matrix for prediction is
#' based on the expected weights of outsample data if they were had been added to
#' insample data to estimate the corresponding MGWRSAR (see Geniaux 2022 for
#' further detail), if method_pred ='shepard'a shepard kernel with k_extra neighbours (default 8) is used and if method_pred='kernel_model' the same kernel
#' and number of neighbors as for computing the MGWRSAR model is used.
#' @return A vector of predictions.
#' @export
#' @rdname predict.mgwrsar
setMethod("predict",'mgwrsar', function(object,newdata, newdata_coords, W = NULL, type = "BPN", h_w = 100,kernel_w = "rectangle",maxobs=4000,beta_proj=FALSE,method_pred='TP', k_extra = 8,...)
{
  z<-predict_mgwrsar(object,newdata, newdata_coords, W , type , h_w,kernel_w ,maxobs,beta_proj,method_pred, k_extra)
  z
}
)

#' Plot method for mgwrsar model
#' @param x   A model of class \code{\link{mgwrsar-class}}.
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
#' @param y   missing
#' @return A Interactive Web Maps with local parameters plot and Open Street Map layer.
#' @aliases plot.mgwrsar
#' @export
#' @method plot mgwrsar
#' @rdname plot.mgwrsar
methods::setMethod("plot", c(x="mgwrsar", y="missing"),function(x, y, type='coef',var=NULL,crs=NULL,mypalette= "RdYlGn",opacity=0.5,fopacity=0.5,nbins=8,radius=500,mytile='Stadia.StamenTonerBackground',myzoom=8,myresolution=150,LayersControl=TRUE,myzoomControl=TRUE,mytile2=NULL,ScaleBar=NULL,ScaleBarOptions=list(maxWidth = 200, metric = TRUE,imperial = FALSE, updateWhenIdle = TRUE),MyLegendTitle=NULL,lopacity=0.5){
                     plot.mgwrsar(x,y,type=type,var=var,crs=crs,mypalette=mypalette,opacity=opacity,fopacity=fopacity,nbins=nbins,radius=radius,mytile=mytile,myzoom=myzoom,myresolution=myresolution,LayersControl=LayersControl,myzoomControl=myzoomControl,mytile2=mytile2,ScaleBar=ScaleBar,ScaleBarOptions=ScaleBarOptions,MyLegendTitle=MyLegendTitle,lopacity=lopacity)
                   })


