#' Class of mgwrsar Model.
#'
#' A S4 class to represent the results of MGWRSAR, TDS-MGWR and related spatial models.
#'
#' @slot Betav matrix. The estimated varying coefficients, dimension (n, kv).
#' @slot Betac numeric. The estimated constant coefficients, length kc.
#' @slot Model character. The type of model (e.g., "GWR", "MGWR", "SAR", "tds_mgwr").
#' @slot fixed_vars character. A vector with names of spatially constant covariates.
#' @slot Y numeric. The dependent variable.
#' @slot XC matrix. The explanatory variables with constant coefficients.
#' @slot XV matrix. The explanatory variables with varying coefficients.
#' @slot X matrix. All explanatory variables.
#' @slot W Matrix. The spatial weight matrix for spatial dependence (row-standardized).
#' @slot isgcv logical. Indicates if Leave-One-Out Cross-Validation (LOOCV) has been computed.
#' @slot edf numeric. The estimated effective degrees of freedom.
#' @slot formula formula. The model formula.
#' @slot data data.frame. The dataframe used for computation.
#' @slot Method character. The estimation technique for spatial dependence models ('2SLS' or 'B2SLS'). Default is '2SLS'.
#' @slot coords matrix. The spatial coordinates of observations.
#' @slot H numeric. The bandwidth vector (spatial).
#' @slot Ht numeric. The bandwidth vector (temporal), if applicable.
#' @slot kernels character. The kernel type(s) used (e.g., 'gauss', 'bisq').
#' @slot adaptive logical. Indicates if an adaptive kernel (nearest neighbors) was used.
#' @slot Type character. The type of Generalized Kernel Product ('GD' for spatial, 'GDT' for spatio-temporal).
#' @slot TP numeric. Indices of target points (if a subset was used).
#' @slot SSRtp numeric. The residual sum of squares calculated only on target points.
#' @slot SSR numeric. The total residual sum of squares.
#' @slot residuals numeric. The vector of residuals.
#' @slot fit numeric. The vector of fitted values.
#' @slot pred numeric. The vector of predicted values (out-of-sample).
#' @slot sev matrix. Local standard errors of varying coefficients.
#' @slot se numeric. Standard errors of constant coefficients.
#' @slot tS numeric. The trace of the Hat matrix (effective number of parameters).
#' @slot Shat matrix. The Hat matrix (or approximation).
#' @slot R_k list. List of partial Hat matrices by covariate (for MGWR inference).
#' @slot h_w numeric. The bandwidth value used for computing the spatial weight matrix W. Default is 0.
#' @slot kernel_w character. The kernel type used for computing W. Default is NULL.
#' @slot RMSE numeric. Root Mean Square Error (on training data).
#' @slot RMSEtp numeric. Root Mean Square Error computed on target points.
#' @slot CV numeric. Leave-One-Out Cross-Validation score.
#' @slot AIC numeric. Akaike Information Criterion.
#' @slot AICc numeric. Corrected Akaike Information Criterion.
#' @slot AICctp numeric. Corrected AIC for target points.
#' @slot BIC numeric. Bayesian Information Criterion.
#' @slot R2 numeric. R-squared.
#' @slot R2_adj numeric. Adjusted R-squared.
#' @slot get_ts logical. Indicates if the trace of the Hat matrix (Tr(S)) was stored.
#' @slot NN numeric. The maximum number of neighbors used for weight computation (truncation parameter).
#' @slot ncore numeric. Number of CPU cores used.
#' @slot mycall call. The original function call.
#' @slot ctime numeric. Computation time in seconds.
#' @slot HRMSE matrix. History of RMSE values (for iterative algorithms like backfitting).
#' @slot HBETA list. History of estimated Beta coefficients at each iteration.
#' @slot loglik numeric. Log-likelihood value.
#' @slot G list. List containing neighboring indices and distances (knn object).
#' @slot V numeric. Sequence of spatial bandwidths tested (for TDS algorithms).
#' @slot Vt numeric. Sequence of temporal bandwidths tested (for TDS algorithms).
#' @slot Z numeric. Temporal or auxiliary variable for GDT kernel type.
#' @slot TS numeric. Diagonal elements of the Hat Matrix.
#' @slot alpha numeric. Ratio parameter for GDT kernels (balancing space and time).
#' @slot HM matrix. Matrix of optimal bandwidths per covariate (for TDS).
#' @slot HKmin numeric. Minimum allowed bandwidth per covariate (for TDS).
#' @slot HKMIN list. List of minimum bandwidths per covariate for spatio-temporal models (TDS).
#' @slot isolated_idx numeric. Indices of observations without sufficient neighbors.
#' @slot my_crs ANY. Coordinate Reference System (CRS) information.
#'
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
           Ht = "numeric",
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
           HM= "matrix",
           HKmin = "numeric",
           HKMIN = "list",
           isolated_idx = "numeric",
           my_crs = "ANY"
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
setMethod("summary", "mgwrsar", function(object, ...) {
  model <- object

  # 1. Verification de la classe
  if (!is(model, 'mgwrsar')) stop("not a mgwrsar object")

  # 2. Extraction des infos
  n <- length(model@Y)
  Model <- model@Model
  Type <- model@Type
  kernels <- model@kernels
  adaptive <- model@adaptive

  # --- EN-TETE ---
  cat("------------------------------------------------------\n")
  cat("Call:\n")
  print(model@mycall)
  cat("\nModel:", Model, "\n")

  # --- KERNELS CONFIGURATION (NOUVEAU BLOC) ---
  cat("------------------------------------------------------\n")
  cat("Kernels Configuration:\n")

  # Fonction pour formater le label (Adaptive vs Fixed)
  get_adapt_label <- function(is_adapt) {
    ifelse(is_adapt, "(Adaptive / Neighbors)", "(Fixed / Distance)")
  }

  # Logique d'affichage selon le Type (Spatial, Temporel ou GDT)
  if (Type == 'GDT') {
    # Spatial
    k_s <- kernels[1]
    a_s <- adaptive[1]
    cat("   Spatial Kernel  :", k_s, get_adapt_label(a_s), "\n")

    # Temporel (Gestion securisee des vecteurs de longueur 1)
    k_t <- if(length(kernels) > 1) kernels[2] else kernels[1]
    a_t <- if(length(adaptive) > 1) adaptive[2] else FALSE # Par defaut Fixed si non specifie
    cat("   Temporal Kernel :", k_t, get_adapt_label(a_t), "\n")

  } else if (Type == 'T') {
    # Temporel pur
    cat("   Temporal Kernel :", kernels[1], get_adapt_label(adaptive[1]), "\n")

  } else {
    # Spatial pur (GD ou defaut)
    cat("   Spatial Kernel  :", kernels[1], get_adapt_label(adaptive[1]), "\n")
  }

  # --- BANDWIDTH CONFIGURATION ---
  cat("------------------------------------------------------\n")
  cat("Bandwidth Configuration:\n")

  simple_models <- c('OLS','GWR','MGWR','MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv',
                     'MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv','MGWRSAR_1_kc_0')
  multiscale_models <- c('tds_mgwr', 'atds_mgwr', 'atds_gwr')

  # Case A: Standard Models (Scalar Bandwidths)
  if (Model %in% simple_models || Model == "OLS") {

    if (Type == 'GD') {
      val <- if(length(model@H) == 1) round(model@H, 2) else "Computed/Optimized"
      cat("   Spatial Bandwidth (H):", val, "\n")

    } else if (Type == 'T') {
      val <- if(length(model@Ht) == 1) round(model@Ht, 2) else "Computed/Optimized"
      cat("   Temporal Bandwidth (H):", val, "\n")

    } else if (Type == 'GDT') {
      h_s <- if(length(model@H) >= 1) round(model@H[1], 2) else "N/A"
      # Handle Ht or second element of H
      h_t_val <-  model@Ht #if(length(model@Ht) == 1) model@Ht else if(length(model@H) > 1) model@H[2] else NA
      h_t <- if(!is.na(h_t_val)) round(h_t_val, 2) else "N/A"

      cat("   Spatial Bandwidth (H) :", h_s, "\n")
      cat("   Temporal Bandwidth (Ht):", h_t, "\n")
    }

    # Case B: Multiscale Models (TDS / ATDS)
  } else if (Model %in% multiscale_models) {

    var_names <- colnames(model@Betav)
    if (is.null(var_names)) var_names <- names(model@H)
    if (is.null(var_names)) var_names <- paste0("Var", 1:length(model@H))

    is_atds <- Model %in% c('atds_mgwr', 'atds_gwr')
    df_bw <- data.frame(Variable = var_names, stringsAsFactors = FALSE)

    # Spatial Bandwidth Column
    df_bw$Spatial_H <- sapply(seq_along(var_names), function(i) {
      if (is_atds && length(model@HM) > 0 && ncol(model@HM) >= i) {
        vals <- na.omit(model@HM[, i])
        if (length(vals) > 1) {
          return(paste0("Adaptive [", round(min(vals), 1), "-", round(max(vals), 1), "]"))
        }
      }
      val <- if(!is.null(names(model@H))) model@H[var_names[i]] else model@H[i]
      if (is.na(val)) "-" else as.character(round(val, 2))
    })

    # Temporal Bandwidth Column
    if (Type == 'GDT') {
      df_bw$Temporal_Ht <- sapply(seq_along(var_names), function(i) {
        val <- NA
        if (length(model@Ht) > 0) {
          val <- if(!is.null(names(model@Ht))) model@Ht[var_names[i]] else model@Ht[i]
        }
        if (length(val) > 1) {
          return(paste0("Adaptive [", round(min(val), 1), "-", round(max(val), 1), "]"))
        } else if (is.na(val)) {
          return("-")
        } else {
          return(as.character(round(val, 2)))
        }
      })
    }

    if (is_atds) cat("   [ATDS] Bandwidth correction during stage 2:\n")
    else cat("   [TDS] Covariate-Specific Bandwidths:\n")
    print(df_bw, row.names = FALSE, right = FALSE)
  }

  # --- MODEL SETTINGS (Suite) ---
  cat("------------------------------------------------------\n")
  cat("Model Settings:\n")
  if (!(Model %in% c('OLS', 'GWR', 'MGWR', 'tds_mgwr', 'atds_mgwr')))
    cat("   Method for spatial autocorrelation:", model@Method, "\n")

  cat("   Computation time:", round(model@ctime, 3), "sec\n")

  tp_check <- if (is.null(model@TP) || length(model@TP) == n) 'NO' else 'YES'
  cat("   Use of Target Points:", tp_check, "\n")
  if (tp_check == 'YES') {
    nb_tp <- if (is.list(model@TP)) length(unlist(model@TP)) else length(model@TP)
    cat("   Number of Target Points:", nb_tp, "\n")
  }

  use_rough <- ifelse(model@NN < n & (!adaptive[1] | (adaptive[1] & model@kernels[1] == 'gauss')),
                      paste0('YES (', model@NN, ' neighbors)'), 'NO')
  cat("   Use of rough kernel approximation:", use_rough, "\n")
  cat("   Parallel computing:", paste0("(ncore = ", model@ncore, ")"), "\n")
  cat("   Number of data points:", n, "\n")

  # --- PARAMETERS ---
  cat("------------------------------------------------------\n")
  cat("Coefficients Summary:\n")

  if (length(model@XC) > 0 || length(model@Betac) > 0) {
    cat("   [Constant Parameters]\n")
    c_names <- names(model@Betac)
    if(is.null(c_names) && length(model@XC) > 0) c_names <- colnames(model@XC)
    vals <- as.numeric(model@Betac)
    names(vals) <- c_names
    print(vals)
    cat("\n")
  }

  if (length(model@XV) > 0) {
    cat("   [Varying Parameters]\n")
    print(summary(model@Betav))
  }

  # --- DIAGNOSTICS ---
  cat("------------------------------------------------------\n")
  cat("Diagnostics:\n")
  if (length(model@edf) > 0) cat("   Effective degrees of freedom:", round(model@edf, 2), "\n")
  if (length(model@tS) > 0 || length(model@AICc) > 0) cat("   AICc:", round(model@AICc, 2), "\n")
  if (length(model@CV) > 0) cat("   LOOCV:", round(model@CV, 4), "\n")

  val_SSR <- if (is.function(model@SSR)) "Not computed" else round(model@SSR, 2)
  val_RMSE <- if (is.function(model@RMSE)) "Not computed" else round(model@RMSE, 4)

  cat("   Residual sum of squares:", val_SSR, "\n")
  cat("   RMSE:", val_RMSE, "\n")
  cat("------------------------------------------------------\n")

  invisible(model)
})
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
#' @param exposant shapenig parameter for tds_mgtwr model, default 6.
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
setMethod("predict",'mgwrsar', function(object,newdata, newdata_coords, W = NULL, type = "BPN", h_w = 100,kernel_w = "rectangle",maxobs=4000,beta_proj=FALSE,method_pred='TP', k_extra = 8,exposant=8,...)
{
  z<-predict_mgwrsar(object,newdata, newdata_coords, W , type , h_w,kernel_w ,maxobs,beta_proj,method_pred, k_extra,exposant)
  z
}
)
#' Plot method for mgwrsar model
#'
#' Maps the results of an mgwrsar model using an interactive Plotly map.
#' Supports plotting coefficients, t-statistics, residuals, and fitted values.
#' Automatically switches between a geographic map (if CRS is present) and a scatter plot.
#'
#' @param x A model of class \code{\link{mgwrsar-class}}.
#' @param y missing (not used).
#' @param type Default 'coef', for plotting the value of the coefficients.
#' Local t-Student can also be plotted using 't_coef', residuals using 'residuals', and fitted values using 'fitted'.
#' @param var Name or index of the variable to plot (required for 'coef' and 't_coef').
#' @param crs A CRS projection (e.g., "EPSG:2154", or an sf/sp object).
#' If NULL, the function attempts to use the CRS stored in the model object.
#' @param mypalette A color palette name (e.g., "RdYlGn", "Viridis", "Portland") or a vector of colors.
#' @param opacity Opacity of the markers (0 to 1).
#' @param size Size of the markers in pixels (default is 5). Replaces the old 'radius' argument.
#' @param title Custom title for the plot. If NULL, a default title is generated.
#' @param show_legend Logical, whether to display the legend/colorbar (default TRUE).
#' @param n_time_steps time step for spatio-temporal model (default 10).
#' @return An interactive Plotly object (Mapbox or Scatter).
#' @aliases plot.mgwrsar
#' @export
#' @importFrom methods setMethod
#' @rdname plot.mgwrsar
methods::setMethod(
  f = "plot",
  signature = c(x = "mgwrsar", y = "missing"),
  definition = function(x, y, type = 'coef', var = NULL, crs = NULL,
                        mypalette = "RdYlGn", opacity = 0.8, size = 5,
                        title = NULL, show_legend = TRUE,n_time_steps = 10,
                        ...) {

    plot.mgwrsar(
      x = x,
      type = type,
      var = var,
      crs = crs,
      mypalette = mypalette,
      size = size,
      opacity = opacity,
      title = title,
      show_legend = show_legend,
      ...
    )
  }
)


