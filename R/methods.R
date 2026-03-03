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
#' @slot edf_k numeric. The estimated effective degrees of freedom of univariate models.
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
#' @slot max_dist numeric. maximum distance possible with this kernel.
#' @slot max_dist_t numeric. maximum temporal distance with this kernel.
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
           edf_k = "numeric",
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
           max_dist="numeric",
           max_dist_t="numeric",
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
#'
#'
#'
#'
setMethod("residuals",'mgwrsar', function(object,...)
{
  object@residuals
}
)
# by_significance_summary_mgwrsar <- function(object,
#                                             alpha_levels = c(0.10, 0.05, 0.01),
#                                             include_fixed = TRUE) {
#
#   get_slot_safe <- function(obj, name) {
#     tryCatch(methods::slot(obj, name), error = function(e) NULL)
#   }
#
#   Betav <- get_slot_safe(object, "Betav")
#   sev   <- get_slot_safe(object, "sev")
#   Betac <- if (include_fixed) get_slot_safe(object, "Betac") else NULL
#   se    <- if (include_fixed) get_slot_safe(object, "se")    else NULL
#   df    <- get_slot_safe(object, "edf")
#   Y     <- get_slot_safe(object, "Y")
#   n     <- length(Y)
#
#   # --- GESTION DU DF MANQUANT (Cas get_AIC = F) ---
#   is_approx_df <- FALSE
#   if (is.null(df) || is.na(df) || length(df) == 0 || df <= 0) {
#     # Approximate with n - 1 if edf is missing, to allow display
#     df <- n - 1
#     is_approx_df <- TRUE
#   }
#
#   has_varying <- !is.null(Betav) && !is.null(sev) && length(Betav) > 0
#   has_fixed   <- !is.null(Betac) && !is.null(se) && length(Betac) > 0
#
#   if (!has_varying && !has_fixed) return(NULL)
#   if (df <= 0) return(NULL)
#
#   row_names <- c("Not significant", "Significant (10% level)",
#                  "Significant ( 5% level)", "Significant ( 1% level)")
#
#   counts_from_pvec <- function(pv) {
#     ok <- is.finite(pv)
#     c(`Not significant` = sum(ok & pv > 0.10),
#       `Significant (10% level)` = sum(ok & pv <= 0.10 & pv > 0.05),
#       `Significant ( 5% level)` = sum(ok & pv <= 0.05 & pv > 0.01),
#       `Significant ( 1% level)` = sum(ok & pv <= 0.01))
#   }
#
#   process_sig <- function(p_matrix) {
#     sum_mat <- matrix(0L, nrow = 4, ncol = ncol(p_matrix),
#                       dimnames = list(row_names, colnames(p_matrix)))
#     for (j in seq_len(ncol(p_matrix))) {
#       sum_mat[, j] <- unname(counts_from_pvec(p_matrix[, j])[row_names])
#     }
#     return(sum_mat)
#   }
#
#   P_RAW <- NULL
#   P_BY  <- NULL
#   nTP <- if(has_varying) nrow(Betav) else 1L
#
#   if (has_varying) {
#     p_raw_v <- 2 * stats::pt(-abs(Betav / sev), df = df)
#     p_by_v <- p_raw_v
#     for (j in seq_len(ncol(p_raw_v))) {
#       p_by_v[, j] <- stats::p.adjust(p_raw_v[, j], method = "BY")
#     }
#     P_RAW <- p_raw_v
#     P_BY  <- p_by_v
#   }
#
#   fixed_counts_raw <- NULL
#   fixed_counts_by  <- NULL
#   if (has_fixed) {
#     p_raw_f <- 2 * stats::pt(-abs(as.numeric(Betac) / as.numeric(se)), df = df)
#     p_by_f <- stats::p.adjust(p_raw_f, method = "BY")
#
#     fixed_counts_raw <- matrix(0L, nrow = 4, ncol = length(Betac), dimnames = list(row_names, names(Betac)))
#     fixed_counts_by  <- matrix(0L, nrow = 4, ncol = length(Betac), dimnames = list(row_names, names(Betac)))
#
#     for(j in seq_along(Betac)) {
#       idx_r <- if(p_raw_f[j] <= 0.01) 4 else if(p_raw_f[j] <= 0.05) 3 else if(p_raw_f[j] <= 0.10) 2 else 1
#       idx_b <- if(p_by_f[j] <= 0.01) 4 else if(p_by_f[j] <= 0.05) 3 else if(p_by_f[j] <= 0.10) 2 else 1
#       fixed_counts_raw[idx_r, j] <- nTP
#       fixed_counts_by[idx_b, j]  <- nTP
#     }
#   }
#
#   res_raw <- if(has_varying) process_sig(P_RAW) else NULL
#   res_by  <- if(has_varying) process_sig(P_BY) else NULL
#
#   if (has_fixed) {
#     res_raw <- if(is.null(res_raw)) fixed_counts_raw else cbind(fixed_counts_raw, res_raw)
#     res_by  <- if(is.null(res_by)) fixed_counts_by else cbind(fixed_counts_by, res_by)
#   }
#
#   return(list(raw = res_raw, adjusted = res_by, is_approx_df = is_approx_df))
# }
by_significance_summary_mgwrsar <- function(object,
                                            alpha_levels = c(0.10, 0.05, 0.01),
                                            include_fixed = TRUE) {

  get_slot_safe <- function(obj, name) {
    tryCatch(methods::slot(obj, name), error = function(e) NULL)
  }

  Betav   <- get_slot_safe(object, "Betav")
  sev     <- get_slot_safe(object, "sev")
  Betac   <- if (include_fixed) get_slot_safe(object, "Betac") else NULL
  se      <- if (include_fixed) get_slot_safe(object, "se")    else NULL

  edf_k_vec <- get_slot_safe(object, "edf_k")
  df_global <- get_slot_safe(object, "edf")
  n         <- length(get_slot_safe(object, "Y"))

  # --- DEFINITION OF IS_APPROX_DF ---
  # Considered as an approximation if the global edf is not available
  is_approx_df <- is.null(df_global) || is.na(df_global) || length(df_global)==0

  has_varying <- !is.null(Betav) && !is.null(sev) && length(Betav) > 0
  has_fixed   <- !is.null(Betac) && !is.null(se) && length(Betac) > 0

  row_names <- c("Not significant", "Significant (10% level)",
                 "Significant ( 5% level)", "Significant ( 1% level)")

  counts_from_pvec <- function(pv) {
    ok <- is.finite(pv)
    c(`Not significant` = sum(ok & pv > 0.10),
      `Significant (10% level)` = sum(ok & pv <= 0.10 & pv > 0.05),
      `Significant ( 5% level)` = sum(ok & pv <= 0.05 & pv > 0.01),
      `Significant ( 1% level)` = sum(ok & pv <= 0.01))
  }

  out <- list(is_approx_df = is_approx_df) # Added here
  nTP <- if(has_varying) nrow(Betav) else 1L

  if (has_varying) {
    p_raw_v <- Betav * 0
    p_by_v  <- Betav * 0

    for (k in colnames(Betav)) {
      df_val <- if(!is.null(edf_k_vec) && !is.na(edf_k_vec[k])) edf_k_vec[k] else
        if(!is.null(df_global) && !is.na(df_global)) df_global else n - 1

      df_val <- max(2, df_val)
      p_raw_v[, k] <- 2 * stats::pt(-abs(Betav[, k] / sev[, k]), df = df_val)
      p_by_v[, k]  <- stats::p.adjust(p_raw_v[, k], method = "BY")
    }

    sum_raw <- matrix(0L, nrow = 4, ncol = ncol(Betav), dimnames = list(row_names, colnames(Betav)))
    sum_by  <- matrix(0L, nrow = 4, ncol = ncol(Betav), dimnames = list(row_names, colnames(Betav)))

    for (k in colnames(Betav)) {
      sum_raw[, k] <- unname(counts_from_pvec(p_raw_v[, k])[row_names])
      sum_by[, k]  <- unname(counts_from_pvec(p_by_v[, k])[row_names])
    }
    out$raw <- sum_raw
    out$adjusted <- sum_by
  }

  if (has_fixed) {
    df_fix <- if(!is.null(df_global) && !is.na(df_global)) df_global else n - 1
    p_raw_f <- 2 * stats::pt(-abs(as.numeric(Betac) / as.numeric(se)), df = max(2, df_fix))
    p_by_f <- stats::p.adjust(p_raw_f, method = "BY")

    fixed_counts_raw <- matrix(0L, nrow = 4, ncol = length(Betac), dimnames = list(row_names, names(Betac)))
    fixed_counts_by  <- matrix(0L, nrow = 4, ncol = length(Betac), dimnames = list(row_names, names(Betac)))

    for(j in seq_along(Betac)) {
      pv_r <- p_raw_f[j]; pv_b <- p_by_f[j]
      idx_r <- if(pv_r <= 0.01) 4 else if(pv_r <= 0.05) 3 else if(pv_r <= 0.10) 2 else 1
      idx_b <- if(pv_b <= 0.01) 4 else if(pv_b <= 0.05) 3 else if(pv_b <= 0.10) 2 else 1
      fixed_counts_raw[idx_r, j] <- nTP
      fixed_counts_by[idx_b, j]  <- nTP
    }
    out$raw <- if(is.null(out$raw)) fixed_counts_raw else cbind(fixed_counts_raw, out$raw)
    out$adjusted <- if(is.null(out$adjusted)) fixed_counts_by else cbind(fixed_counts_by, out$adjusted)
    out$is_approx_df<-is_approx_df
  }

  return(out)
}
#' summary for mgwrsar model
#'
#' @param object A model of class \code{\link{mgwrsar-class}}.
#' @param ... summary parameters forwarded.
#' @return A summary object.
#' @export
#' @rdname summary.mgwrsar
setMethod("summary", "mgwrsar", function(object, ...) {
  model <- object
  if (!is(model, 'mgwrsar')) stop("not a mgwrsar object")

  n <- length(model@Y)
  Model <- model@Model
  Type <- model@Type
  kernels <- model@kernels
  adaptive <- model@adaptive

  cat("------------------------------------------------------\n")
  cat("Call:\n")
  print(model@mycall)
  cat("\nModel:", Model, "\n")

  cat("------------------------------------------------------\n")
  cat("Kernels Configuration:\n")
  get_adapt_label <- function(is_adapt) ifelse(is_adapt, "(Adaptive / Neighbors)", "(Fixed / Distance)")

  if (Type == 'GDT') {
    cat("    Spatial Kernel  :", kernels[1], get_adapt_label(adaptive[1]), "\n")
    k_t <- if(length(kernels) > 1) kernels[2] else kernels[1]
    a_t <- if(length(adaptive) > 1) adaptive[2] else FALSE
    cat("    Temporal Kernel :", k_t, get_adapt_label(a_t), "\n")
  } else if (Type == 'T') {
    cat("    Temporal Kernel :", kernels[1], get_adapt_label(adaptive[1]), "\n")
  } else {
    cat("    Spatial Kernel  :", kernels[1], get_adapt_label(adaptive[1]), "\n")
  }

  cat("------------------------------------------------------\n")
  cat("Bandwidth Configuration:\n")
  simple_models <- c('OLS','GWR','MGWR','MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv','MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv','MGWRSAR_1_kc_0')
  multiscale_models <- c('tds_mgwr', 'atds_mgwr', 'atds_gwr')

  if (Model %in% simple_models || Model == "OLS") {
    if (Type == 'GD') {
      cat("    Spatial Bandwidth (H):", if(length(model@H) == 1) round(model@H, 2) else "Optimized", "\n")
    } else if (Type == 'T') {
      cat("    Temporal Bandwidth (H):", if(length(model@Ht) == 1) round(model@Ht, 2) else "Optimized", "\n")
    } else if (Type == 'GDT') {
      cat("    Spatial Bandwidth (H) :", if(length(model@H) >= 1) round(model@H[1], 2) else "N/A", "\n")
      cat("    Temporal Bandwidth (Ht):", if(!is.null(model@Ht)) round(model@Ht, 2) else "N/A", "\n")
    }
  } else if (Model %in% multiscale_models) {
    var_names <- if(!is.null(colnames(model@Betav))) colnames(model@Betav) else names(model@H)
    df_bw <- data.frame(Variable = var_names, Spatial_H = sapply(var_names, function(vn) round(model@H[vn], 2)))
    if (Type == 'GDT') df_bw$Temporal_Ht <- sapply(var_names, function(vn) round(model@Ht[vn], 2))
    cat(paste0("    [", toupper(Model), "] Covariate-Specific Bandwidths:\n"))
    print(df_bw, row.names = FALSE)
  }

  cat("------------------------------------------------------\n")
  cat("Model Settings:\n")
  cat("    Computation time:", round(model@ctime, 3), "sec\n")
  cat("    Number of data points:", n, "\n")

  cat("------------------------------------------------------\n")
  cat("Coefficients Summary:\n")
  if (length(model@Betac) > 0) {
    cat("    [Constant Parameters]\n")
    print(model@Betac)
    cat("\n")
  }
  if (length(model@Betav) > 0) {
    cat("    [Varying Parameters]\n")
    print(summary(model@Betav))
  }

  # --- DUAL SIGNIFICANCE ANALYSIS (raw and BY-adjusted) ---
  if (length(model@sev) > 0) {
    sig_res <- by_significance_summary_mgwrsar(model)
    if (!is.null(sig_res)) {
      if (isTRUE(sig_res$is_approx_df)) {
        cat("\n    [Significance Analysis - WARNING: Local Approximation]\n")
        cat("    Note: edf not computed (get_AIC=F). Using df = n-1 for local p-values.\n")
      } else {
        cat("\n    [Significance Analysis - Raw (Unadjusted)]\n")
      }
      print(sig_res$raw)
      cat("\n    [Significance Analysis - BY Adjusted (FDR Control)]\n")
      print(sig_res$adjusted)
      cat("\n    Note: Counts reflect number of points with significant local coefficients.\n")
    }
  }

  cat("------------------------------------------------------\n")
  cat("Diagnostics:\n")
  if (length(model@edf) > 0 && !is.na(model@edf)) {
    cat("    Effective degrees of freedom:", round(model@edf, 2), "\n")
  } else {
    cat("    Effective degrees of freedom: Not computed (get_AIC=F)\n")
  }
  if (length(model@AICc) > 0 && !is.na(model@AICc)) cat("    AICc:", round(model@AICc, 2), "\n")
  cat("    RMSE:", if (is.function(model@RMSE)) "Not computed" else round(model@RMSE, 4), "\n")
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


