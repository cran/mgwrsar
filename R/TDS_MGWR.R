#' Top-Down Scale (TDS) and Adaptive Top-Down Scale (ATDS) Estimation for MGWR
#'
#' @description
#' This function implements the "Top-Down Scale" (TDS) methodology for estimating
#' Multiscale Geographically Weighted Regression (MGWR) models.
#' Unlike classical backfitting approaches that fully optimize bandwidths at each iteration, TDS uses a pre-defined sequence of decreasing bandwidths
#' to efficiently identify the optimal spatial scale for each covariate.
#'
#' The function supports two main algorithms:
#' \itemize{
#'   \item \strong{'tds_mgwr'}: A backfitting algorithm that selects a unique optimal bandwidth for each covariate from a decreasing sequence.
#'   \item \strong{'atds_mgwr'}: Extends 'tds_mgwr' with a second "boosting" stage (Adaptive TDS). It refines estimates by allowing bandwidths to vary locally, capturing complex spatial patterns (e.g., simultaneous large-scale trends and local variations).
#' }
#'
#' @usage TDS_MGWR(formula, data, coords, Model = 'tds_mgwr',
#'                     kernels = 'gauss', fixed_vars = NULL, Ht = NULL,
#'                     control_tds = list(nns = 25, get_AIC = FALSE, init_model = "OLS"),
#'                     control = list(adaptive = TRUE))
#'
#' @param formula A formula object specifying the model (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing the variables in the model.
#' @param coords A matrix or data frame of coordinates (2 columns for spatial).
#' @param Model A character string specifying the algorithm. Options:
#' \itemize{
#'   \item \code{'tds_mgwr'} (default): Top-Down Scale MGWR (Stage 1 only).
#'   \item \code{'atds_mgwr'}: Adaptive Top-Down Scale MGWR (Stage 1 + Stage 2 boosting).
#'   \item \code{'atds_gwr'}: Univariate Adaptive Top-Down Scale GWR.
#' }
#' @param kernels A character string or vector specifying the kernel type (e.g., \code{'triangle'}, \code{'bisq'}, \code{'gauss'}). Default is \code{'triangle'}.
#' @param fixed_vars A character vector indicating the names of variables with spatially stationary (fixed) coefficients. Default is \code{NULL}.
#' @param Ht Numeric. Optional bandwidth for the second dimension (time) if using spatio-temporal models (Type 'GDT').
#' @param control_tds A named list of control parameters specific to the TDS algorithm:
#' \describe{
#'   \item{\code{nns}}{Integer. Number of bandwidth steps in the decreasing sequence (default 30).}
#'   \item{\code{get_AIC}}{Logical. If \code{TRUE}, computes AICc (slower). Default \code{FALSE} (except for \code{'atds_mgwr'}).}
#'   \item{\code{init_model}}{Character. Initial model type to start backfitting: \code{'OLS'} (default), \code{'GWR'}, \code{'GTWR'}, or \code{'known'}.}
#'   \item{\code{ncore}}{Integer. Number of cores for parallelization. Default 1.}
#'   \item{\code{tol}}{Numeric. Convergence tolerance. Default 0.001.}
#'   \item{\code{nrounds}}{Integer. Number of boosting rounds for Stage 2 (only for \code{'atds_mgwr'}). Default 3.}
#' }
#' @param control A named list of standard control arguments passed to the internal \code{MGWRSAR} calls:
#' \describe{
#'   \item{\code{adaptive}}{Logical or Vector. \code{TRUE} for adaptive bandwidth (nearest neighbors), \code{FALSE} for fixed distance. Can be a vector of length 2 for space/time.}
#'   \item{\code{Type}}{Character. Spatial weighting type: \code{'GD'} (Spatial only) or \code{'GDT'} (Space-Time).}
#'   \item{\code{NN}}{Integer. Maximum number of neighbors for distance matrix computation (truncation). Default is \code{nrow(data)}.}
#' }
#'
#' @details
#' The TDS algorithm works in two stages:
#' \enumerate{
#'   \item **Stage 1 (Backfitting):** Starts with a global model (OLS) or a simple GWR. It iteratively updates the bandwidth for each covariate by testing values from a decreasing sequence. This avoids the "yo-yo" effect of standard backfitting and converges faster.
#'   \item **Stage 2 (Boosting - atds_mgwr only):** Uses the residuals from Stage 1 to iteratively refine coefficients. This stage allows the effective bandwidth to adapt locally, improving accuracy for covariates with spatially heterogeneous scales.
#' }
#'
#' @return An object of class \code{mgwrsar} containing:
#' \item{Betav}{Matrix of spatially varying coefficients.}
#' \item{H}{Vector of optimal bandwidths found for each covariate.}
#' \item{RMSE}{Root Mean Square Error of the final model.}
#' \item{AICc}{Corrected Akaike Information Criterion (if requested).}
#' \item{fitted.values}{Vector of fitted values.}
#' \item{residuals}{Vector of residuals.}
#'
#' @references
#' Geniaux, G. (2024). Top-Down Scale Approaches for Multiscale GWR with Locally Adaptive Bandwidths. \emph{Springer Nature}.
#'
#' @seealso \code{\link{MGWRSAR}}, \code{\link{golden_search_bandwidth}}
#' @export
TDS_MGWR <- function(formula, data, coords,
                         Model = 'tds_mgwr',
                         kernels = 'gauss',
                         fixed_vars = NULL,
                         Ht = NULL,
                         control_tds = list(nns = 25, get_AIC = FALSE, init_model = "OLS"),
                         control = list(adaptive = TRUE)) {

  set_bandwidth_bounds <- function(coords, X, time = NULL,
                                   adaptive = c(FALSE, FALSE),
                                   quant_spatial = c(0.05, 0.95),
                                   quant_temporal = c(0.05, 0.95),
                                   scale_time = FALSE,
                                   scale_factor = NULL) {
    # ----------------------------------------------------------
    # Checks and formatting
    # ----------------------------------------------------------
    if (!is.matrix(coords) && !is.data.frame(coords))
      stop("coords must be a matrix or data.frame with 2 columns (x, y)")
    coords <- as.matrix(coords)
    if (ncol(coords) < 2)
      stop("coords must have at least two columns (x, y)")

    n <- nrow(coords)
    K <- ncol(X)

    # Harmonize adaptive parameter
    if (length(adaptive) == 1)
      adaptive <- rep(adaptive, 2)
    names(adaptive) <- c("spatial", "temporal")

    message("-> Setting bandwidth ranges (adaptive = ",
            paste(adaptive, collapse = ", "), ")")

    # ----------------------------------------------------------
    # [1] Spatial range
    # ----------------------------------------------------------
    if (isTRUE(adaptive["spatial"])) {
      # adaptive = TRUE -> bounds in number of neighbors
      hs_range <- c(K + 2, n)
    } else {
      # non-adaptive -> bounds based on distances
      nsample <- min(2000, n)
      idx <- sample(seq_len(n), nsample)
      D_sp <- as.matrix(dist(coords[idx, ]))
      D_sp <- D_sp[lower.tri(D_sp)]
      hs_min <- quantile(D_sp, quant_spatial[1], na.rm = TRUE)
      hs_max <- quantile(D_sp, quant_spatial[2], na.rm = TRUE)
      hs_range <- c(hs_min, hs_max)
    }

    # ----------------------------------------------------------
    # [2] Temporal range (if applicable)
    # ----------------------------------------------------------
    if (!is.null(time)) {
      if (isTRUE(adaptive["temporal"])) {
        ht_range <- c(K + 2, n)
      } else {
        if (any(is.na(time))) stop("`time` must not contain NA values")
        time <- as.numeric(time)
        D_t <- dist(time)
        D_t <- as.numeric(D_t)
        ht_min <- quantile(D_t, quant_temporal[1], na.rm = TRUE)
        ht_max <- quantile(D_t, quant_temporal[2], na.rm = TRUE)
        ht_range <- c(ht_min, ht_max)

        if (isTRUE(scale_time)) {
          if (is.null(scale_factor)) {
            nsample <- min(2000, n)
            idx <- sample(seq_len(n), nsample)
            D_sp <- as.numeric(dist(coords[idx, ]))
            ratio_st <- median(D_sp, na.rm = TRUE) / median(D_t, na.rm = TRUE)
            scale_factor <- ratio_st
            message(sprintf("    Scaling time dimension by factor = %.3f", scale_factor))
          }
          ht_range <- ht_range * scale_factor
        }
      }
    } else {
      ht_range <- NULL
    }

    # ----------------------------------------------------------
    # [3] Conditional output (GWR vs GTWR)
    # ----------------------------------------------------------
    if (is.null(time) && length(adaptive) == 1) {
      # GWR case only -> no temporal dimension
      out <- list(hs_range = hs_range)
    } else {
      out <- list(hs_range = hs_range, ht_range = ht_range)
      if (exists("scale_factor")) out$scale_factor <- scale_factor
    }

    # ----------------------------------------------------------
    # [4] Console summary
    # ----------------------------------------------------------
    message(sprintf("    hs_range = [%.3f, %.3f]", hs_range[1], hs_range[2]))
    if (!is.null(ht_range))
      message(sprintf("    ht_range = [%.3f, %.3f]", ht_range[1], ht_range[2]))

    return(out)
  }

  # ============================================================
  #    TDS-MGWR Main Routine
  #    Clean, modular, and unified version
  # ============================================================
  if(is.null(control$Type)) control$Type='GD'
  message("\n-------------------------------------------")
  if(control$Type=='GD')  message("Running Top-Down Scale MGWR") else if(control$Type=='GDT')  message("Running Top-Down Scale MGTWR")
  message("-------------------------------------------")
  start <- proc.time()

  # ============================================================
  # 1. INITIALIZATION
  # ============================================================
  init_param_tds <- function(env = parent.frame()) {
    with(env, {

      `%||%` <- function(x, y) if (!is.null(x)) x else y

      mf <- model.frame(formula, data)
      data <- data[, names(mf)]

      # ============================================================
      # [1] MODEL CHECKS AND BASIC DIMENSIONS
      # ============================================================
      if (!(Model %in% c('tds_mgwr', 'atds_mgwr', 'atds_gwr')))
        stop('Only atds_gwr, tds_mgwr and atds_mgwr Model can be estimated using Top Down Scale approach in this release.')

      n_time <- m <- n <- nrow(data)

      # ============================================================
      # [2] INITIALIZATION OF TDS PARAMETERS
      # ============================================================
      if (is.null(control_tds$randomize_order)) control_tds$randomize_order=TRUE
      type <- 'proportional'
      if (is.null(control$Type)) control$Type <- 'GD'
      if (is.null(control_tds$init_model) & control$Type == 'GD') control_tds$init_model <- 'OLS'
      if (is.null(control_tds$init_model) & control$Type == 'GDT') control_tds$init_model <- 'GTWR'

      if (is.null(control_tds$nns)) control_tds$nns <- 30
      if (is.null(control_tds$blocksize)) control_tds$blocksize <- n
      if (is.null(control_tds$get_AIC)) control_tds$get_AIC <- FALSE
      if (Model == 'atds_mgwr') control_tds$get_AIC <- TRUE
      if (is.null(control_tds$ncore)) control_tds$ncore <- 1
      if (is.null(control_tds$tol)) control_tds$tol <- 0.001
      if (is.null(control_tds$extra_iter)) control_tds$extra_iter <- 3
      if (is.null(control_tds$refine)) control_tds$refine <- FALSE
      if (is.null(control_tds$check_pairs)) control_tds$check_pairs <- FALSE
      if (is.null(control_tds$maxit)) control_tds$maxit <- 100


      if (is.null(control_tds$nrounds)) control_tds$nrounds <- 3
      if (is.null(control_tds$verbose)) control_tds$verbose <- FALSE

      if (is.null(control_tds$V)) V <- control_tds$V <- NULL
      if (is.null(control_tds$min_dist)) min_dist <- NULL
      if (is.null(control_tds$first_nn)) control_tds$first_nn <- n
      if (is.null(control_tds$BETA)) control_tds$BETA <- NULL
      if (is.null(control_tds$TRUEBETA)) control_tds$TRUEBETA <- NULL
      if (is.null(control_tds$H)) H <- control_tds$H <- NULL
      if (is.null(control_tds$browser)) control_tds$browser <- 0

      # ============================================================
      # [3] PARALLELIZATION SAFETY AND CONTROL PARAMETERS
      # ============================================================
      if (is.null(control$NN)) control$NN <- n
      if (is.null(control$isgcv)) control$isgcv <- FALSE
      if (is.null(control$Type)) control$Type <- 'GD'
      if (is.null(control$TP)) control$TP <- 1:n
      if (is.null(control$alpha)) control$alpha <- 1

      control$get_ts <- TRUE
      control$get_s <- control_tds$get_AIC

      # ============================================================
      # [4] COORDINATE HANDLING AND DUPLICATE CHECK
      # ============================================================
      if(control$Type %in% c('GD','GDT')) coords<-make_unique_by_structure(coords)
      if(control$Type %in% c('GDT','T'))  control$Z<-make_unique_by_structure(control$Z)
      # ============================================================
      # [5] DISTANCE MATRICES PREPARATION
      # ============================================================
      if (!is.null(control$Z)) coords_in <- as.matrix(cbind(coords, control$Z)) else coords_in <- coords
      if (!('indexG' %in% names(control)) & is.null(control_tds$model_stage1)) {
        G <- prep_d(coords = coords_in, NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
        control$indexG <- G$indexG
        control$dists <- G$dists
      } else if (!is.null(model_init)) {
        if (is.null(model_init@G))
          G <- prep_d(coords = coords_in, NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
        else
          G <- model_init@G
        control$indexG <- G$indexG
        control$dists <- G$dists
      } else {
        G <- list()
        G$indexG <- control$indexG
        G$dists <- control$dists
      }
      if (control$Type=='GDT') {
        controlvd <- modifyList(control, list(dists = NULL, indexG = NULL, Type = "GD",adaptive = control$adaptive[1]))
        Gtemp <- prep_d(coords = coords_in[,1:2], NN = control$NN, TP = control$TP, kernels = kernels[1], Type = 'GD')
        controlvd$indexG <- Gtemp$indexG
        controlvd$dists <- Gtemp$dists

        controlvt <- modifyList(control, list( Type = "T",dists = NULL, indexG = NULL,adaptive = FALSE))
        Gtemp <- prep_d(coords = as.matrix(coords_in[,3],ncol=1), NN = control$NN, TP = control$TP, kernels = kernels[2], Type = 'T')
        controlvt$indexG <- Gtemp$indexG
        controlvt$dists <- Gtemp$dists
        rm(Gtemp)
      }


      if (control$adaptive[1]) max_dist <- n
      else max_dist <- max(control$dists[[1]][, ncol(control$dists[[1]])])

      # ============================================================
      # [6] MODEL FRAME PREPARATION
      # ============================================================
      mycall <- match.call()
      mf <- model.frame(formula, data)
      mt <- attr(x = mf, which = "terms")
      Y <- model.extract(mf, "response")
      X <- model.matrix(object = mt, data = mf)
      colnames(X) <- clean_colnames(X)

      if (!is.null(fixed_vars)) fixed_vars <- intersect(fixed_vars, colnames(X))
      if (ncol(X) > 1)
        "tds_gwr is only designed for univariate regression. With multivariate GWR, you cannot be sure of achieving a global optimum. The results may show better in-sample RMSE or AICc compared to GWR, tds_mgwr, or atds_mgwr, but it must be compared to these models using cross-validation."

      formula <- formula(lm(formula, data))

      # ============================================================
      # [7] PREVIOUS MODEL HANDLING (STAGE 1)
      # ============================================================
      if (!is.null(control_tds$model_stage1)) {
        S <- control_tds$model_stage1@Shat
        if (!(control_tds$model_stage1@Model %in% c('OLS', 'GWR')))
          V <- control_tds$V <- control_tds$model_stage1@V
        BETA <- control_tds$model_stage1@Betav
        myAICc <- control_tds$model_stage1@AICc
      }

      # ============================================================
      # [8] MATRICES AND DESIGN INFO
      # ============================================================
      I <- Diagonal(n)
      namesX <- colnames(X)
      if (colnames(X)[1] == "(Intercept)") colnames(X)[1] <- namesX[1] <- 'Intercept'
      data$Intercept <- 1
      K <- length(namesX)

      if (is.null(control_tds$minv)) {
        if (kernels[1] == 'gauss') control_tds$minv <- minv <- 2 else if (is.null(control_tds$minv)) control_tds$minv <- minv <- K + 1
      }
      if (!('BETA' %in% ls())) BETA <- NULL
      if (!('TRUEBETA' %in% ls())) TRUEBETA <- NULL
      varying <- setdiff(namesX, fixed_vars)

      # ============================================================
      # [9] FORMULAE FOR LOCAL AND CONSTANT PARTS
      # ============================================================
      myformula_b <- as.formula(paste0('e0~-1+', paste0(varying, collapse = '+')))
      if (length(fixed_vars) > 0)
        formula_constant <- as.formula(paste0('e0~-1+', paste0(fixed_vars, collapse = '+')))

      # ============================================================
      # [10] ALGORITHM HISTORY INITIALIZATION
      # ============================================================
      HBETA <- list()
      HBETA[[1]] <- BETA
      if (control_tds$get_AIC) {
        HAICc <- c()
        HTS <- list()
      }
      HOPT <- rep(NA, K)
      names(HOPT) <- namesX
    })
  }
  init_param_tds()
  reassign_control(control_tds)

  TSik <- matrix(0, nrow = n, ncol = length(varying))
  # Remove TSik updates in update_opt etc.
  colnames(TSik) = namesX
  idx_init <- idx <- 1:n


  # ============================================================
  # 2. BUILD BANDWIDTH SEQUENCES
  # ============================================================
  built_Vseq <- function(env = parent.frame()) {
    with(env, {

      # ============================================================
      # [0] Initialization and defaults
      # ============================================================
      temporal_distance_modulo <- function(x, cycling = 365) {
        x_mod <- x %% cycling
        x_mod[x_mod == 0] <- cycling
        pmin(x_mod, cycling - x_mod)
      }

      if (is.null(first_nn)) first_nn <- n
      if (is.null(control$NN)) control$NN <- n
      first_nn <- min(first_nn, control$NN)

      if (is.null(control$adaptive)) control$adaptive <- c(TRUE, FALSE)
      if (length(control$adaptive) == 1) control$adaptive <- rep(control$adaptive, 2)

      # ============================================================
      # [1] Build spatial sequence V
      # ============================================================
      if (is.null(V)) {
        if (type == "proportional") {
          alpha2 <- exp(log(minv / first_nn) / (nns + 1))
          V <- sapply(1:(nns + 1), function(x) round(first_nn * (alpha2)^x))
        } else {
          v <- round(first_nn / nns)
          V <- seq(first_nn, minv, by = -v)
        }
        V <- unique(V[V >= minv])
      } else {
        m <- V[1]
        V <- V[V < (n - 2)]
      }

      # ============================================================
      # [2] Filter by minimum distance or NN thresholds
      # ============================================================
      if (!is.null(minv)) V <- V[V > minv]

      if (Model != "atds_gwr") V <- c(first_nn, V)

      if (is.null(control_tds$min_dist) && !control$adaptive[1]) {
        min_dist <- min(V)
      } else if (!is.null(control_tds$min_dist) && !control$adaptive[1]) {
        V <- V[V > control_tds$min_dist]
        min_dist <- control_tds$min_dist
      }

      l <- length(V)
      V5 <- V[unique(round(quantile(seq_along(V), c(1, 0.75, 0.5, 0.25, 0))))]

      # ============================================================
      # [3] Temporal candidate sequence Vt (if GDT and non-adaptive)
      # ============================================================
      Vt <- NULL
      min_dist_t <- NULL
      max_dist_t <- NULL

      if (isTRUE(control$Type == "GDT") && !control$adaptive[2]) {

        kernels_t <- unlist(strsplit(kernels[2], "_"))[1]
        format_t  <- unlist(strsplit(kernels[2], "_"))[2]
        cycling   <- as.numeric(unlist(strsplit(kernels[2], "_"))[3])

        # Select correct temporal distance matrix
        if (!is.na(cycling) && "dist_t_modulo" %in% names(G$dists)) {
          dist_t <- abs(G$dists[["dist_t_modulo"]])
        } else {
          dist_t <- abs(G$dists[["dist_t"]])
        }

        if (!is.na(cycling)) {
          max_dist_t <- round(cycling / 2)
        } else {
          max_dist_t <- max(dist_t, na.rm = TRUE)
        }

        min_dist_t <- quantile(dist_t[dist_t > 0.0000001], 0.02, na.rm = TRUE)
        alpha2 <- exp(log(min_dist_t / max_dist_t) / (nns + 1))
        Vt <- sapply(1:(nns + 1), function(x) round(max_dist_t * (alpha2)^x))
        Vt <- Vt[Vt >= min_dist_t]
        Vt <- unique(c(max_dist_t, Vt))
        V5t <- Vt[unique(round(quantile(seq_along(Vt), c(1, 0.75, 0.5, 0.25, 0))))]

        if (is.null(control_tds$min_dist_t)) {
          min_dist_t <- min(Vt)
        } else {
          Vt <- Vt[Vt > control_tds$min_dist_t]
          V5t <- V5t[V5t > control_tds$min_dist_t]
          min_dist_t <- control_tds$min_dist_t
        }
      }

      # ============================================================
      # [4] Spatial distances and conversion according to adaptivity
      # ============================================================
      if (isTRUE(control$adaptive[1])) {
        max_dist <- min(n, control$NN)
        min_dist <- min(V)
      } else {
        # Non-adaptive kernels use true distances
        if (!"dist_s" %in% names(G$dists))
          stop("Missing 'dist_s' in G$dists for spatial distances")

        max_dist <- max(G$dists[["dist_s"]], na.rm = TRUE)
        # Convert neighbour indices (V) into actual distances
        V <- c(max_dist, sapply(V, function(x) median(G$dists[["dist_s"]][, x], na.rm = TRUE)))
        min_dist <- min(V)
      }
    })
  }
  built_Vseq()


  # ============================================================
  # 2b. GET V Vt pairs with sufficient local variance
  # ============================================================

  get_HKmin <- function(env = parent.frame()) {
    cat(
      "\n------------------------------------------------------------\n",
      "Finding minimum bandwidths by covariate\n",
      "Criterion: quantile(local_var / global_var) > threshold\n",
      "------------------------------------------------------------\n",
      sep = ""
    )

    res <- with(env, {

      ## 1. Global settings ----
      min_var_ratio <- if (!is.null(control$min_var_ratio)) control$min_var_ratio else 0.05
      min_var_q     <- if (!is.null(control$min_var_q))     control$min_var_q     else 0.10

      # indices of "varying" covariates in X
      varying_idx <- match(varying, colnames(X))
      if (any(is.na(varying_idx))) {
        stop("Some names in 'varying' are not found in 'colnames(X)'.")
      }

      # Pre-calculation of global variance per covariate (for all varying)
      var_glob_all <- vapply(
        varying_idx,
        function(j) stats::var(X[, j]),
        FUN.VALUE = numeric(1)
      )

      near_zero <- (!is.finite(var_glob_all) | var_glob_all <= .Machine$double.eps)

      HKmin <- NULL   # for Type = "GD"
      HKMIN <- NULL   # for Type = "GDT"

      ## 2. Type = "GD": pure spatial bandwidths ----
      if (identical(control$Type, "GD")) {

        V_sorted <- sort(unique(V))

        # HKmin: 1 value per covariate in 'varying'
        HKmin <- rep(NA_real_, length(varying))
        names(HKmin) <- varying

        # we do not treat the intercept (position 1) in the variance test,
        # we will fill it later as min(HKmin[-1])
        cov_idx_to_check <- which(!near_zero)   # among all varying
        cov_idx_to_check <- setdiff(cov_idx_to_check, 1L)  # exclude Intercept

        l <- 1L
        while (length(cov_idx_to_check) > 0L && l <= length(V_sorted)) {

          v <- V_sorted[l]

          # spatial weights for this v
          stage1 <- prep_w(
            H        = c(v),
            kernels = kernels,
            Type     = control$Type,
            adaptive = control$adaptive,
            dists    = G$dists,
            indexG   = G$indexG,
            alpha    = 1
          )
          Wd <- stage1$W

          # Sub-matrix of X for "varying" covariates
          X_var <- X[, varying_idx, drop = FALSE]

          # Local means and local variances for all varying at once
          Mu   <- as.matrix(Wd %*% X_var)             # n x p_varying
          Mu2  <- as.matrix(Wd %*% (X_var^2))         # n x p_varying
          VarL <- Mu2 - Mu^2                          # n x p_varying

          # local/global ratios per variable (columns)
          # VarL[ , k] / var_glob_all[k]
          ratio_mat <- sweep(VarL, 2L, var_glob_all, FUN = "/")

          # quantile of ratios per variable
          q_ratio_all <- apply(
            ratio_mat,
            2L,
            stats::quantile,
            probs    = min_var_q,
            na.rm    = TRUE,
            names    = FALSE
          )

          # among covariates not yet fixed (excluding intercept),
          # check which ones pass the threshold for this v
          for (kpos in cov_idx_to_check) {
            q_r <- q_ratio_all[kpos]
            if (is.finite(q_r) && q_r > min_var_ratio) {
              HKmin[kpos] <- v
            }
          }

          # Update list of covariates to check
          cov_idx_to_check <- cov_idx_to_check[is.na(HKmin[cov_idx_to_check])]

          l <- l + 1L
        }

        # Intercept: assign min of HKmin (excluding NA and intercept)
        if (length(HKmin) > 1L) {
          HKmin[1L] <- min(HKmin[-1L], na.rm = TRUE)
        }

        # fallback for covariates that never reached the threshold
        if (any(is.na(HKmin))) {
          warning("Some variables never reached the minimal variance ratio; ",
                  "setting HKmin to max(V) for those variables.")
          HKmin[is.na(HKmin)] <- max(V_sorted)
        }

        if (any(HKmin >= max_dist)) {
          cat(names(HKmin[HKmin >= max_dist]), "\n")
          cat(HKmin[HKmin >= max_dist], "\n")
          stop("remove covariates with HKmin == max_dist)")
        }
      }

      ## 3. Type = "GDT": spatio-temporal bandwidths (hs, ht) ----
      if (identical(control$Type, "GDT")) {

        hs_vals <- sort(unique(V))
        ht_vals <- sort(unique(Vt))
        n_hs <- length(hs_vals)
        n_ht <- length(ht_vals)

        p_vary <- length(varying_idx)

        HKMIN <- vector("list", p_vary)
        names(HKMIN) <- varying

        valid_vars <- which(!near_zero)
        valid_vars <- setdiff(valid_vars, 1L)

        if (length(valid_vars) > 0L) {

          X_var         <- X[, varying_idx[valid_vars], drop = FALSE]
          var_glob_sub <- var_glob_all[valid_vars]

          # --- NEW: warm start on ht ---
          # j_start = 1 at the beginning (smallest ht)
          j_start <- 1L

          for (ii in seq_len(n_hs)) {

            hs <- hs_vals[ii]

            # for this hs, the first acceptable ht
            ht_for_var <- rep(NA_real_, length(valid_vars))

            # ht loop starting at the smallest known necessary ht
            for (jj in j_start:n_ht) {       # <--- main acceleration here

              ht <- ht_vals[jj]

              stage1 <- prep_w(
                H        = c(hs, ht),
                kernels = kernels,
                Type     = control$Type,
                adaptive = control$adaptive,
                dists    = G$dists,
                indexG   = G$indexG,
                alpha    = 1
              )
              Wd <- stage1$W

              Mu_sub    <- as.matrix(Wd %*% X_var)
              Mu2_sub  <- as.matrix(Wd %*% (X_var^2))
              VarL_sub <- Mu2_sub - Mu_sub^2

              ratio_sub <- sweep(VarL_sub, 2L, var_glob_sub, "/")
              q_ratio_sub <- apply(ratio_sub, 2L, quantile,
                                   probs = min_var_q, na.rm = TRUE, names = FALSE)

              not_yet <- which(is.na(ht_for_var))

              if (length(not_yet) > 0L) {
                newly_ok <- not_yet[
                  is.finite(q_ratio_sub[not_yet]) &
                    (q_ratio_sub[not_yet] > min_var_ratio)
                ]

                if (length(newly_ok) > 0L) {
                  ht_for_var[newly_ok] <- ht
                }
              }

              if (all(!is.na(ht_for_var))) break
            }

            # determine the minimal ht retained for this hs
            if (any(!is.na(ht_for_var))) {
              j_start <- min(match(ht_for_var, ht_vals), na.rm = TRUE)
            }

            # record the boundary
            for (idx_loc in which(!is.na(ht_for_var))) {
              kpos_global <- valid_vars[idx_loc]
              HKMIN[[kpos_global]] <- rbind(
                HKMIN[[kpos_global]],
                c(hs = hs, ht = ht_for_var[idx_loc])
              )
            }
          }
        }

        if (p_vary >= 2L && !is.null(HKMIN[[2L]])) {
          HKMIN[[1L]] <- HKMIN[[2L]]
        }
      }
      list(HKmin = HKmin, HKMIN = HKMIN)
    })

    # inject back into calling environment
    if (!is.null(res$HKmin)) {
      assign("HKmin", res$HKmin, envir = env)
    }
    if (!is.null(res$HKMIN)) {
      assign("HKMIN", res$HKMIN, envir = env)
    }

    invisible(NULL)
  }

  if(control_tds$check_pairs) get_HKmin()


  # ============================================================
  # 3. STARTING MODEL
  # ============================================================
  starting_model <- function(env = parent.frame()) {
    with(env, {
      message("Initializing starting model...")

      # Case 1: Pre-computed model (provided in control_tds$model_known)
      if (!is.null(control_tds$init_model) && control_tds$init_model == "known") {
        message("-> Using pre-computed model from control_tds$model_known")

        if (is.null(control_tds$model_known))
          stop("control_tds$model_known must be provided when init_model='known'")

        mod0 <- control_tds$model_known
        BETA0 <- mod0$BETA
        residuals0 <- mod0$residuals

        if (isTRUE(control_tds$get_AIC)) {
          S0 <- mod0$S %||% NULL
          Rk0 <- mod0$Rk %||% NULL
        } else {
          S0 <- Rk0 <- NULL
        }
      } else {

        # Case 2: Model to estimate
        init_type <- control_tds$init_model
        message(" Estimating starting model of type: ", init_type)

        if (init_type == "OLS") {
          message("    Fitting global OLS model...")
          model0 <- MGWRSAR(formula = formula, data = data, coords = coords, fixed_vars = NULL, Model = 'OLS', H = NULL, kernels = NULL, control = list())
          BETA0 = matrix(coef(model0)$Betac, byrow = TRUE, nrow = nrow(data), ncol = K)
          colnames(BETA0) = namesX
          residuals0 <- model0@residuals
          if (control$adaptive[1]) model0@H = n else model0@H = max_dist
          H <- rep(n, length(namesX))
          names(H) = namesX
          if (control$Type == 'GDT') {
            Ht <- model0@Ht <- max_dist_t
          }

          if (control_tds$get_AIC) {
            Rkk <<- Rk <- list()
            XXtX <- eigenMapMatMult(solve(crossprod(X)), t(X))
            rownames(XXtX) <- colnames(X)
            S =  eigenMapMatMult(X, XXtX)
            for (k in namesX) {
              Rk[[k]] <- outer(X[, k], XXtX[k, ], '*')
            }
            model0@TS <- diag(S)
            myAICc = model0@AIC
          }
        }
        if (init_type %in% c("GWR")) {
          message("    Fitting spatial GWR/MGWR model...")
          controlv <- control
          controlv$adaptive <- controlv$adaptive[1]
          controlv$Type = 'GD'
          myrange <- set_bandwidth_bounds(coords, X, time = control$Z, quant_spatial = c(0.005, 0.999), adaptive = control$adaptive)
          lower.bound <- myrange$hs_range[1]
          upper.bound <- myrange$hs_range[2]
          if (control_tds$get_AIC) {
            controlv$get_s = TRUE
            controlv$get_Rk = TRUE
          }
          res <- golden_search_bandwidth(formula = formula, Ht = NULL, data = data, coords = coords, fixed_vars = fixed_vars, kernels = kernels[1], Model = 'GWR', control = controlv, lower.bound = lower.bound, upper.bound = upper.bound)
          cat(" optimal bandwidth = ", res$minimum)
          model0 = res$model
          BETA0 = model0@Betav
          H <- rep(model0@H, length(namesX))
          names(H) = namesX
          myAICc = model0@AIC
          S <- model0@Shat
          Rk <- model0@R_k ## NULL
          if (control$Type == 'GDT') {
            Ht <- model0@Ht <- max_dist_t
          }
        }

        if (init_type %in% c("GTWR")) {
          message("    Fitting spatio-temporal GTWR/MGTWR model...")

          myrange <- set_bandwidth_bounds(coords, X, time = control$Z, quant_spatial = c(0.005, 0.999), adaptive = control$adaptive)
          control_tds_temp = control_tds
          if (control_tds$get_AIC) control_tds_temp$get_AIC = TRUE

          res_st <- search_bandwidths(
            formula = formula,
            data = data,
            coords = coords,
            kernels = kernels,#,
            Model = "GWR",
            control = control,
            hs_range = myrange$hs_range,
            ht_range = myrange$ht_range,
            n_seq = 10,
            ncore = control_tds$ncore,
            n_rounds = 1
          )
          model0 = res_st$best_model
          BETA0 = model0@Betav
          H <- rep(model0@H, length(namesX))
          names(H) = namesX
          if (control$Type == 'GDT') {
            Ht <- model0@Ht
          }
          myAICc = model0@AIC
          S <- model0@Shat
          Rk <- model0@R_k
        }

        if (init_type %in% c("MGWR")) {
          message("    Fitting spatio-temporal GTWR/MGTWR model...")
          control_tds_temp <- control_tds
          control_tds_temp$init_model = 'OLS'
          if (control_tds$get_AIC) control_tds_temp$get_AIC = TRUE
          model0 <- TDS_MGWR(formula, data, coords, control_tds_temp, control)
          S <- model0@Shat
          Rk <- model0@R_k
          BETA0 = model0@Betav
          myAICc = model0@AIC
          H = model0@H
          if (control$Type == 'GDT') {
            Ht <- model0@Ht <- max_dist_t
          }
        }

        BETA = BETA0
        data$e0 <- e0 <- model0@residuals

      }

    })
  }

  starting_model()

  # ============================================================
  # 4. STAGE 1 - MAIN BACKFITTING
  # ============================================================

  init_bandwidth_bounds_from_model <- function(env = parent.frame()) {
    with(env, {
      # ------------------------------------------------------------
      # [1]  INITIALIZE DEFAULTS
      # ------------------------------------------------------------
      if (control$Type == 'GDT') {
        ddown  <- down  <- up  <- opt  <- model0@H
        if(control_tds$init_model != 'GTWR') last_opt_t <- ddown_t <- down_t <- up_t <- opt_t <- rep(max_dist_t, lvarying) else last_opt_t <- ddown_t <- down_t <- up_t <- opt_t <- rep(model_lm0@Ht, lvarying)
        if (control_tds$test)
          opt_t <- rep(control_tds$tH, lvarying)
      } else {
        up        <- rep(V[max(i - 1, 1)], lvarying)
        last_opt <- opt <- rep(V[i], lvarying)
        ddown     <- down <- rep(V[min(i + 1, length(V))], lvarying)
        opt_t     <- up_t <- down_t <- ddown_t <- NULL
      }

      # ------------------------------------------------------------
      # [2]  ADAPT INITIAL VALUES ACCORDING TO STARTING MODEL
      # ------------------------------------------------------------
      if (!is.null(model0)) {
        model_type <- model0@Model

        if (model_type == "GWR") {
          opt[1:K] <- model0@H
          if(control$Type == 'GDT')  opt_t[1:K] <- max_dist_t
        } else if (model_type == "MGWR") {
          opt <- model0@H
          if(control$Type == 'GDT')  opt_t[1:K] <- max_dist_t
        } else if (model_type == "MGTWR") {
          opt    <- model0@H
          opt_t <- model0@Ht
        } else if(model_type == "OLS"){
          if(control$adaptive[1]) opt[1:K] <- min(control$NN,n) else opt[1:K] <- max_dist
          if(control$Type == 'GDT')   opt_t[1:K] <- max_dist_t
        }
      }

      # ------------------------------------------------------------
      # [3]  COMPUTE SPATIAL BOUNDS
      # ------------------------------------------------------------
      for (k in seq_len(K)) {

        # Upper bound
        if (opt[k] < max_dist)
          up[k] <- tail(V[V > opt[k]], 1)
        else
          up[k] <- max_dist

        # Lower and double-lower bounds
        if (opt[k] > min_dist) {
          if (sum(V < opt[k]) > 0)
            down[k] <- head(V[V < opt[k]], 1)
          else
            down[k] <- min_dist

          ddown[k] <- max(
            min(down[down > min_dist], ddown[k], na.rm = TRUE),
            min_dist, na.rm = TRUE
          )
        } else {
          ddown[k] <- down[k] <- min_dist
        }
      }

      # ------------------------------------------------------------
      # [4]  COMPUTE TEMPORAL BOUNDS (IF APPLICABLE)
      # ------------------------------------------------------------
      if (control$Type == 'GDT' && !is.null(opt_t)) {
        for (k in seq_len(K)) {

          # Upper bound
          if (opt_t[k] < max_dist_t)
            up_t[k] <- tail(Vt[Vt > opt_t[k]], 1)
          else
            up_t[k] <- max_dist_t

          # Lower and double-lower bounds
          if (opt_t[k] > min_dist_t) {
            if (sum(Vt < opt_t[k]) > 0)
              down_t[k] <- head(Vt[Vt < opt_t[k]], 1)
            else
              down_t[k] <- min_dist_t

            ddown_t[k] <- max(
              min(down_t[down_t > min_dist_t], ddown_t[k], na.rm = TRUE),
              min_dist_t, na.rm = TRUE
            )
          } else {
            ddown_t[k] <- down_t[k] <- min_dist_t
          }
        }
      }
      # Label bounds if relevant
      names(stable) <- names(up) <- names(opt) <- names(down) <- names(ddown) <- varying
      if (control$Type == 'GDT')
        names(up_t) <- names(opt_t) <- names(down_t) <- names(ddown_t) <- varying
    })
  }

  stage1_tds_mgwr_old <- function(env = parent.frame()) {
    with(env, {

      # ============================================================
      #  INITIALIZATION & PRE-PROCESSING
      # ============================================================
      OPT = TRUE

      # Blocksize default for large datasets
      if (is.null(control_tds$blocksize) & n >= 4000)
        control_tds$blocksize = 500

      # Create target points partition
      TP <- quadTP(coords, control_tds$blocksize)
      foldsl <- split(seq_len(nrow(TP)), TP$id)

      # ------------------------------------------------------------
      # 1.1 Initialize algorithm iteration index and starting bandwidths
      # ------------------------------------------------------------
      i = 1
      if (init_model == 'GWR') {
        ts <- model_lm0@tS
      } else {
        new_ts <- ts <- model_lm0@tS
        new_TS <- model_lm0@TS
      }


      # ------------------------------------------------------------
      # 1.2 Verbose startup message
      # ------------------------------------------------------------
      if (verbose) {
        if (!is.null(model_stage1)) {
          cat('Starting from a previous model : \n')
          summary(model_stage1)
        } else {
          cat('\n i =', i, ' Starting model \n')
        }
      }

      # ------------------------------------------------------------
      # 1.3 Initialize basic metrics and control variables
      # ------------------------------------------------------------
      lvarying = length(varying)
      rmse <- sqrt(mean(data$e0^2))
      stable <- rep(0, length(varying))
      spacestable = TRUE
      switched = FALSE



      if (verbose)
        cat(paste0(
          '\n\n########################################################################\n',
          ' STAGE 1 : find unique bandwidth for each covariate  \n',
          ' using Top Down Scale Approach with backfitting for Type = ', control$Type,
          '  ...\n########################################################################\n'
        ))

      # ============================================================
      # 2 INITIAL BANDWIDTH CONFIGURATION
      # ============================================================
      init_bandwidth_bounds_from_model()

      if(any(!is.na(control_tds$H))) {
        for(i in 1:lvarying) {
          opt[varying[i]] <- control_tds$H[i]
        }
        stable = stable + 8
      }
      if(any(!is.na(control_tds$Ht))) {
        for(i in 1:lvarying) {
          opt_t[varying[i]] <- control_tds$Ht[i]
        }
        stable = stable + 8
      }

      # ------------------------------------------------------------
      # 2.1 Initialize AIC-related tracking structures (if needed)
      # ------------------------------------------------------------
      if (control_tds$get_AIC) {
        BestCrit <- lastCrit <- lastCrit2 <- 10^6
        AIC_deep <- NA
        Last_best_AICg <- last_AICc <- myAICc + 10^6
      }

      # ------------------------------------------------------------
      # 2.2 Parallelization setup
      # ------------------------------------------------------------
      if (Sys.info()[['sysname']] == "Linux") {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
        try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
      }

      if (ncore)
         registerDoParallel(cores = ncore)
      else
         registerDoSEQ()

      op <- if (control_tds$ncore>1) "%dopar%" else "%do%"
      fop <- get(op, envir = asNamespace("foreach"))

      # Initialize best parameters and convergence deltas
      bestBETA = BETA
      delta_rmse = 1
      delta_AICc = 1

      # ============================================================
      # [3] BACKFITTING MAIN LOOP
      # ============================================================
      if (control_tds$test) nnrep = 2 else nnrep = 5
      extra_iter <- control_tds$extra_iter
      post_conv_count <- 0
      bestRMSE = rmse
      converged = FALSE
      while ((post_conv_count < extra_iter || i < control_tds$test ) & i<control_tds$maxit) {

        # ------------------------------------------------------------
        # 3.1 Iteration bookkeeping
        # ------------------------------------------------------------
        last_rmseG = rmse
        BETA = bestBETA
        if (verbose) cat('\n\n##############\n ', i)


        if (control_tds$get_AIC) {
          last_AICc = myAICc
          BestCrit = lastCrit
          St = S
        }

        varyingT <- varying
        if (i == browser) browser()
        gc()

        myformula_b = as.formula(paste0('e0~-1+', paste0(varyingT, collapse = '+')))
        controlv <- control
        last_opt = opt
        if (control$Type == 'GDT') last_opt_t = opt_t

        # ------------------------------------------------------------
        # 3.2 MAIN LOOP OVER VARYING VARIABLES
        # ------------------------------------------------------------
        for (k in varying) {

          if (verbose) cat(' ', k)
          last_rmse = rmse

          # --- Bandwidth update depending on model type ---
          if(!converged){
            if (any(stable < ifelse(control_tds$refine, 2, 4))) {
              if (control$Type == 'GD') update_opt()
              if (control$Type == 'GDT' & !control_tds$test) update_opt_st()
              if (control$Type == 'GDT' & control_tds$test) update_opt_st_test()
            } else if(control_tds$refine) {
              if (control$Type == 'GD'  ) {
                cat(
                  '\n one step Golden search ratio \n')
                if(k == tail(varying, 1)) converged = TRUE
                data$e0k <- data$e0 + BETA[, k] * X[, k]
                myformula_bk = as.formula(paste0('e0k~-1+', k))

                refined <- golden_search_bandwidth( formula = myformula_bk, Ht = NULL, data = data, coords = coords,
                                                    fixed_vars = NULL, kernels = kernels, Model = "GWR", control = control,
                                                    lower.bound = HKmin[k], upper.bound = V[max(1, which(V == opt[k]) - 3)])
                model_k <- refined$model
                opt[k] <- model_k@H
                betav <- model_k@Betav
                e0 <- residuals(model_k)
                isol <- is.na(betav)
                e0[isol] <- data$e0[isol]
                if (get_AIC) {
                  TSik[!isol, k] <- model_k@TS[!isol]
                  Sk <- model_k@Shat
                }
              }
              # else if (control$Type == 'GDT') {
              #   if(k==tail(varying,1)) converged=TRUE
              #   data$e0k<-data$e0+BETA[,k]*X[,k]
              #   myformula_bk=as.formula(paste0('e0k~-1+',k))
              #
              #   refined <- golden_search_2d_bandwidth(formula = myformula_bk,
              #                                         data = data, coords = coords, fixed_vars = NULL,
              #                                         kernels = kernels, Model =  'GWR', control = control,
              #                                         lower.bound.space = min(V), upper.bound.space = V[max(1,which(V==opt[k])-3)],
              #                                         lower.bound.time = min(Vt), upper.bound.time = Vt[max(1,which(Vt==opt_t[k])-3)], tolerance_s =1, tolerance_t=1)
              #
              #   model_k<-refined$model
              #   opt[k]<-model_k@H
              #   opt_t[k]<-model_k@Ht
              #   betav<-model_k@Betav
              #   e0 <- residuals(model_k)
              #   isol <- is.na(betav)
              #   e0[isol] <- data$e0[isol]
              #   if (get_AIC) {
              #     TSik[!isol, k] <- model_k@TS[!isol]
              #     Sk <- model_k@Shat
              #   }
              # }
            } else  update_opt_known()
          } else            {
            update_opt_known()
          }

          # --- Update AIC trace and matrices ---            # --- Update local coefficients ---

          if (control_tds$get_AIC) {
            Rkk[[k]] <- compute_Rk(Rk[[k]], Sk, St, foldsl)
            St = St + Rkk[[k]] - Rk[[k]]
            new_TS = diag(St)
            new_ts = sum(new_TS)

            if (verbose & control_tds$get_AIC)
              cat(' AICc : ', aicc_f(e0[idx_init], new_ts, n_time))

            # if(  aicc_f(e0[idx_init], new_ts, n_time)<Last_best_AICg) {
            #   BETA[!isol, k] = betav[!isol]
            #   Last_best_AICg<-aicc_f(e0[idx_init], new_ts, n_time)
            # }

          }
          BETA[!isol, k] = betav[!isol]

          # --- Update residuals and RMSE ---
          fit = rowSums(BETA * X)
          data$e0 = Y - fit
          rmse = sqrt(mean(data$e0^2))
        }
        # ------------------------------------------------------------
        # 3.3 FIXED VARIABLES UPDATE
        # ------------------------------------------------------------
        for (k in fixed_vars) {
          last_rmse = rmse
          data$e0k <- data$e0 + BETA[, k] * X[, k]
          myformula_bk = as.formula(paste0('e0k~-1+', k))

          model_tds <- MGWRSAR(formula = myformula_bk, data = data,
                               coords = coords, fixed_vars = k, kernels = kernels,
                               H = NULL, Model = 'OLS', control = controlv)

          BETA[, k] = model_tds@Betac

          if (control_tds$get_AIC) {
            TSik[, k] <- unlist(res[mybest, 'TS'])
            Sk <- unlist(res[mybest, 'S'][[1]])
            Rkk[[k]] <- compute_Rk(Rk[[k]], Sk, St, foldsl)
            St = St + Rkk[[k]] - Rk[[k]]
            new_TS = diag(St)
            new_ts = sum(new_TS)
          }

          fit = rowSums(BETA * X)
          data$e0 = Y - fit
          rmse = sqrt(mean(data$e0^2))
        }

        # ------------------------------------------------------------
        # 3.4 UPDATE CRITERIA AND CONVERGENCE CHECKS
        # ------------------------------------------------------------
        if (control_tds$get_AIC)
          myAICc <- aicc_f(data$e0[idx_init], new_ts, n_time)

        # Update stability tracker
        if (control$Type == 'GDT') {
          stable[opt != last_opt | opt_t != last_opt_t] <- 0
          stable[opt == last_opt & opt_t == last_opt_t] <- stable[opt == last_opt & opt_t == last_opt_t] + 1
        } else {
          stable[opt != last_opt] <- 0
          stable[opt == last_opt] <- stable[opt == last_opt] + 1
        }

        # Compute deltas
        delta_rmse = (last_rmseG - sqrt(mean((Y - rowSums(BETA * X))^2))) / last_rmseG
        if (control_tds$get_AIC)
          delta_AICc = (last_AICc - myAICc) / last_AICc

        # ------------------------------------------------------------
        # 3.5 IF IMPROVEMENT: SAVE BEST STATE
        # ------------------------------------------------------------
        if (control_tds$get_AIC) {
          S = St
          for (k in varying)
            Rk[[k]] <- Rkk[[k]]
        }

        if (sqrt(mean((Y - rowSums(BETA * X))^2)) <= bestRMSE) {
          bestBETA = BETA
          bestRMSE = sqrt(mean((Y - rowSums(BETA * X))^2))
          H = opt
          if (control$Type == 'GDT') Ht = opt_t
          mybestG = i + 1
          if (control_tds$get_AIC) {
            mybestS = S
            mybestRk = Rk
          }
        }

        # ------------------------------------------------------------
        # 3.6 END-OF-ITERATION HOUSEKEEPING
        # ------------------------------------------------------------
        fit = rowSums(BETA * X)
        data$e0 = Y - fit
        rmse = sqrt(mean(data$e0^2))
        delta_rmse = (last_rmseG - rmse) / last_rmseG

        if (control_tds$get_AIC) {
          HTS[[i + 1]] <- new_TS
          HAICc <- c(HAICc, myAICc)
        }

        HBETA[[i + 1]] <- BETA

        if (!is.null(TRUEBETA)) {
          for (k in 1:K)
            HRMSE[i + 1, k] = sqrt(mean((TRUEBETA[, k] - BETA[, k])^2))
          HRMSE[i + 1, K + 1] <- mean(HRMSE[i + 1, 1:K])
          HRMSE[i + 1, K + 2] <- opt[k]
          if (control_tds$get_AIC) HRMSE[i + 1, K + 3] <- myAICc
          HRMSE[i + 1, K + 4] <- sqrt(mean(data$e0^2))
        }

        if (verbose & control_tds$get_AIC)
          cat('\n delta_AICc: ', delta_AICc, ' AICc : ', myAICc, ' stable ', stable, '\n delta_rmse ', delta_rmse, ' rmse ', rmse)
        if (verbose & !control_tds$get_AIC)
          cat('\n stable ', stable, ' delta_rmse ', delta_rmse, ' rmse ', rmse)
        if (verbose) cat('\n\n H =', opt)
        if (verbose & control$Type == 'GDT') cat('\n Ht =', opt_t)


        # --- sign_history to detect oscillations ---
        if (!exists("sign_history")) sign_history <- integer(0)
        if (!exists("delta_history")) delta_history <- numeric(0)
        if (!exists("bandwidth_history")) bandwidth_history <- list()
        if (!exists("pp_streak")) pp_streak <- 0L  # compteur ping-pong parfait

        current_sign <- sign(delta_rmse)
        if (!is.na(current_sign) && current_sign != 0) {
          sign_history <- c(sign_history, current_sign)
          delta_history <- c(delta_history, delta_rmse)
        }

        if (control$Type == "GDT")
          bandwidth_history[[length(bandwidth_history) + 1]] <- list(h = opt, ht = opt_t)
        else
          bandwidth_history[[length(bandwidth_history) + 1]] <- list(h = opt)

        # Limit sign_history size
        if (length(sign_history) > 10) {
          sign_history <- tail(sign_history, 10)
          delta_history <- tail(delta_history, 10)
        }
        if (length(bandwidth_history) > 10)
          bandwidth_history <- tail(bandwidth_history, 10)

        # ============================================================
        # [1] Detect strict Yoyo  (RMSE + - + - +)
        # ============================================================
        if (length(sign_history) >= 4) {
          diffs <- diff(sign_history)
          # (+ - + -)
          alt_seq <- all(abs(diffs) == 2) && all(diff(diffs) != 0)
          alt_count <- if (alt_seq) length(diffs[diffs != 0]) else 0
        } else {
          alt_count <- 0
        }

        # Mean of tow last RMSE
        mean_last2 <- if (length(delta_history) >= 2)
          mean(abs(tail(delta_history, 2))) else abs(delta_rmse)

        # ============================================================
        # [2]  Detect ping-pong between bandwidths
        # ============================================================
        bw_pingpong <- 0
        if (length(bandwidth_history) >= 4) {
          for (j in seq_len(length(bandwidth_history) - 2)) {
            bw_prev <- bandwidth_history[[j + 1]]$h
            bw1 <- bandwidth_history[[j]]$h
            bw2 <- bandwidth_history[[j + 2]]$h

            has_change <- !isTRUE(all.equal(bw_prev, bw1, tolerance = 0))

            ## Checks that it returns to the original value two iterations later
            same_h <- isTRUE(all.equal(bw1, bw2, tolerance = 0))

            same_ht <- TRUE
            if (control$Type == "GDT") {
              bw_prev_t <- bandwidth_history[[j + 1]]$ht
              bw1_t <- bandwidth_history[[j]]$ht
              bw2_t <- bandwidth_history[[j + 2]]$ht

              has_change <- has_change || !isTRUE(all.equal(bw_prev_t, bw1_t, tolerance = 0))
              same_ht <- isTRUE(all.equal(bw1_t, bw2_t, tolerance = 0))
            }

            # ping-pong count
            if (has_change && same_h && same_ht)
              bw_pingpong <- bw_pingpong + 1
          }
        }

        # ============================================================
        # [3] 3 times of perfect ping-pong
        # ============================================================
        detect_perfect_pingpong <- function(h1, h2, h3) {
          isTRUE(all.equal(h3, h1, tolerance = 0)) && !isTRUE(all.equal(h3, h2, tolerance = 0))
        }

        if (length(bandwidth_history) >= 3) {
          bw_tm2 <- bandwidth_history[[length(bandwidth_history) - 2L]]
          bw_tm1 <- bandwidth_history[[length(bandwidth_history) - 1L]]
          bw_t   <- bandwidth_history[[length(bandwidth_history)]]

          pp_spatial <- detect_perfect_pingpong(bw_tm2$h, bw_tm1$h, bw_t$h)
          pp_temporal <- TRUE
          if (control$Type == "GDT")
            pp_temporal <- detect_perfect_pingpong(bw_tm2$ht, bw_tm1$ht, bw_t$ht)

          if (pp_spatial && pp_temporal)
            pp_streak <- pp_streak + 1L
          else
            pp_streak <- 0L

          if (pp_streak >= 3L) {
            cat("\n[WARNING] Perfect ping-pong detected 3 times in a row -- inserting intermediate bandwidths.\n")

            # ----- SPATIAL -----
            h_now  <- as.numeric(bw_t$h)
            h_prev <- as.numeric(bw_tm1$h)

            if (isTRUE(control$adaptive[1])) {
              need_mid_h <- abs(h_now - h_prev) > 2
              h_mid <- round((h_now + h_prev) / 2)
            } else {
              rel_diff_h <- abs(h_now - h_prev) / pmax(1e-12, (abs(h_now) + abs(h_prev)) / 2)
              need_mid_h <- rel_diff_h > 0.001
              h_mid <- (h_now + h_prev) / 2
            }

            for (kk in seq_along(h_now)) {
              if (isTRUE(need_mid_h[kk])) {
                V <- sort(unique(c(V, h_mid[kk])))
                opt[kk] <- h_mid[kk]
              }
            }

            # ----- TEMPORAL -----
            if (control$Type == "GDT") {
              ht_now  <- as.numeric(bw_t$ht)
              ht_prev <- as.numeric(bw_tm1$ht)

              if (isTRUE(control$adaptive[2])) {
                need_mid_ht <- abs(ht_now - ht_prev) > 2
                ht_mid <- round((ht_now + ht_prev) / 2)
              } else {
                rel_diff_ht <- abs(ht_now - ht_prev) / pmax(1e-12, (abs(ht_now) + abs(ht_prev)) / 2)
                need_mid_ht <- rel_diff_ht > 0.001
                ht_mid <- (ht_now + ht_prev) / 2
              }

              for (kk in seq_along(ht_now)) {
                if (isTRUE(need_mid_ht[kk])) {
                  Vt <- sort(unique(c(Vt, ht_mid[kk])))
                  opt_t[kk] <- ht_mid[kk]
                }
              }
            }

            cat("# Intermediate bandwidths inserted. Resuming iterations.\n")
            pp_streak <- 0L
            sign_history <- integer(0)
            delta_history <- numeric(0)
            bandwidth_history <- list(list(h = opt, ht = if (control$Type == "GDT") opt_t else NULL))
          }
        }

        # ============================================================
        # [4] stopping conditions
        # ============================================================
        stop_flag <- FALSE
        reason <- NULL

        # (a) weak yoyo oscillation
        if (alt_count >= 3 && mean_last2 < tol) {
          stop_flag <- TRUE
          reason <- sprintf("weak yoyo oscillation (alt=%d, mean|RMSE|=%.6f < tol=%.6f)",
                            alt_count, mean_last2, tol)
        }

        # (b) repetition stricte des memes bandwidths
        if (bw_pingpong >= 3) {
          stop_flag <- TRUE
          reason <- sprintf("repeated bandwidth ping-pong (%d identical alternations)", bw_pingpong)
        }

        if (stop_flag) {
          cat("\n Stopping criterion reached:", reason, "\n")

          opt <- round((opt + last_opt) / 2)
          if (control$Type == "GDT")
            opt_t <- round((opt_t + last_opt_t) / 2)
          i = i + 1
          break
        }

        if (abs(delta_rmse) < tol) post_conv_count <- post_conv_count + 1 else post_conv_count <- 0

        i = i + 1

      }

      # ============================================================
      # [4] FINALIZATION & MODEL ASSEMBLY
      # ============================================================
      names(H) <- varying

      if (control_tds$get_AIC) {
        myAICc <- HAICc[mybestG - 1]
        TS = as.numeric(unlist(HTS[mybestG]))
        tS = sum(TS)
      }

      BETA <- bestBETA
      if (!is.null(TRUEBETA))
        HRMSE <- HRMSE[1:mybestG, ]
      HBETA <- HBETA[1:mybestG]

      fit = rowSums(BETA * X)
      data$e0 = Y - fit

      # ------------------------------------------------------------
      # 4.1 BUILD FINAL MODEL OBJECT
      # ------------------------------------------------------------
      modelGWR <- new('mgwrsar')
      modelGWR@mycall <- mycall
      modelGWR@data <- data
      modelGWR@coords <- as.matrix(coords)
      modelGWR@X <- modelGWR@XV <- X
      if (verbose) cat('\n Time = ', (proc.time() - start)[3], '\n')

      modelGWR@formula = formula
      modelGWR@Model = Model
      modelGWR@H = H
      modelGWR@fixed_vars <- unique(c(as.character(fixed_vars),
                                      names(modelGWR@H)[which(modelGWR@H == max_dist)]))
      modelGWR@Betav = BETA
      modelGWR@fit = rowSums(BETA * X)
      modelGWR@residuals = Y - modelGWR@fit
      modelGWR@RMSE = sqrt(mean(modelGWR@residuals^2))
      modelGWR@Type = control$Type
      modelGWR@kernels = kernels
      modelGWR@adaptive = control$adaptive
      modelGWR@TP = 1:n
      modelGWR@ncore = control_tds$ncore
      modelGWR@NN = control$NN

      if(control_tds$check_pairs) {
        if(!is.null(HKmin))  modelGWR@HKmin <- HKmin
        if(!is.null(HKMIN)) modelGWR@HKMIN <- HKMIN
      }
      if (!is.null(control$Z))
        modelGWR@Z = control$Z
      modelGWR@V = V

      if (control$Type == 'GDT') {
        modelGWR@Vt = Vt
        modelGWR@Ht = Ht
        modelGWR@alpha = control$alpha
      }

      modelGWR@Y = Y
      modelGWR@ctime <- (proc.time() - start)[3]
      modelGWR@HBETA <- HBETA

      if (control_tds$get_AIC) {
        if (!is.null(S)) modelGWR@Shat = mybestS
        modelGWR@TS = diag(mybestS)
        modelGWR@tS = sum(modelGWR@TS)
        modelGWR@edf <- n - tS
        modelGWR@AICc <- myAICc
        for (k in varying)
          modelGWR@R_k[[k]] <- mybestRk[[k]]
      }

      if (!is.null(TRUEBETA))
        modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[, 1]), ]

      if(Model == 'atds_mgwr') modelGWR@G = G
      control_tds$model_stage1 <- model_stage1 <- returned_model <- modelGWR
    })
  }

  stage1_tds_mgwr <- function(env = parent.frame()) {
    with(env, {

      # ============================================================
      #  INITIALIZATION & PRE-PROCESSING
      # ============================================================
      OPT = TRUE

      # Blocksize default for large datasets
      if (is.null(control_tds$blocksize) & n >= 4000)
        control_tds$blocksize = 500

      # Create target points partition
      TP <- quadTP(coords, control_tds$blocksize)
      foldsl <- split(seq_len(nrow(TP)), TP$id)

      # ------------------------------------------------------------
      # 1.1 Initialize algorithm iteration index and starting bandwidths
      # ------------------------------------------------------------
      i = 1
      if (init_model == 'GWR') {
        ts <- model_lm0@tS
      } else {
        new_ts <- ts <- model_lm0@tS
        new_TS <- model_lm0@TS
      }

      # ------------------------------------------------------------
      # 1.2 Verbose startup message
      # ------------------------------------------------------------
      if (verbose) {
        if (!is.null(model_stage1)) {
          cat('Starting from a previous model : \n')
          summary(model_stage1)
        } else {
          cat('\n i =', i, ' Starting model \n')
        }
      }

      # ------------------------------------------------------------
      # 1.3 Initialize basic metrics and control variables
      # ------------------------------------------------------------
      lvarying = length(varying)
      rmse <- sqrt(mean(data$e0^2))
      stable <- rep(0, length(varying))
      spacestable = TRUE
      switched = FALSE

      # Definition of the detection function (defined once here)
      detect_perfect_pingpong <- function(h1, h2, h3) {
        isTRUE(all.equal(h3, h1, tolerance = 0)) && !isTRUE(all.equal(h3, h2, tolerance = 0))
      }

      if (verbose)
        cat(paste0(
          '\n\n########################################################################\n',
          ' STAGE 1 : find unique bandwidth for each covariate  \n',
          ' using Top Down Scale Approach with backfitting for Type = ', control$Type,
          '  ...\n########################################################################\n'
        ))

      # ============================================================
      # 2 INITIAL BANDWIDTH CONFIGURATION
      # ============================================================
      init_bandwidth_bounds_from_model()

      if(any(!is.na(control_tds$H))) {
        for(i in 1:lvarying) {
          opt[varying[i]] <- control_tds$H[i]
        }
        stable = stable + 8
      }
      if(any(!is.na(control_tds$Ht))) {
        for(i in 1:lvarying) {
          opt_t[varying[i]] <- control_tds$Ht[i]
        }
        stable = stable + 8
      }

      # ------------------------------------------------------------
      # 2.1 Initialize AIC-related tracking structures
      # ------------------------------------------------------------
      if (control_tds$get_AIC) {
        BestCrit <- lastCrit <- lastCrit2 <- 10^6
        AIC_deep <- NA
        Last_best_AICg <- last_AICc <- myAICc + 10^6
      }

      # ------------------------------------------------------------
      # 2.2 Parallelization setup
      # ------------------------------------------------------------
      if (Sys.info()[['sysname']] == "Linux") {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
        try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
      }

      if (ncore>1)
        registerDoParallel(cores = ncore)
      else
        registerDoSEQ()

      op <- if (control_tds$ncore>1) "%dopar%" else "%do%"
      fop <- get(op, envir = asNamespace("foreach"))

      # Initialize best parameters and convergence deltas
      bestBETA = BETA
      delta_rmse = 1
      delta_AICc = 1

      # Tracking variables for oscillations
      if (!exists("sign_history")) sign_history <- integer(0)
      if (!exists("delta_history")) delta_history <- numeric(0)
      if (!exists("bandwidth_history")) bandwidth_history <- list()
      if (!exists("pp_streak")) pp_streak <- 0L

      # ============================================================
      # [3] BACKFITTING MAIN LOOP
      # ============================================================
      if (control_tds$test) nnrep = 5 else nnrep = 5
      extra_iter <- control_tds$extra_iter
      post_conv_count <- 0
      bestRMSE = rmse
      converged = FALSE
      patience=0
      while ((post_conv_count < extra_iter || i < control_tds$test ) & i < control_tds$maxit) {

        # ------------------------------------------------------------
        # 3.1 Iteration bookkeeping
        # ------------------------------------------------------------
        last_rmseG = rmse
        BETA = bestBETA
        if (verbose) cat('\n\n ', i)

        if (control_tds$get_AIC) {
          last_AICc = myAICc
          BestCrit = lastCrit
          St = S
        }

        varyingT <- varying
        if (i == browser) browser()
        gc()

        myformula_b = as.formula(paste0('e0~-1+', paste0(varyingT, collapse = '+')))
        controlv <- control
        last_opt = opt
        if (control$Type == 'GDT') last_opt_t = opt_t

        # ------------------------------------------------------------
        # 3.2 MAIN LOOP OVER VARYING VARIABLES
        # ------------------------------------------------------------
        if(control_tds$randomize_order) {
          set.seed(as.integer(Sys.time()))
          samp_varying=varying[c(1,sample(2:length(varying),length(varying)-1))]
        } else samp_varying

        for (k in samp_varying) {

          if (verbose) cat(' ', k)
          last_rmse = rmse

          # --- Bandwidth update depending on model type ---
          if(!converged){
            if (any(stable < ifelse(control_tds$refine, 2, 4))) {
              if (control$Type == 'GD') update_opt()
              if (control$Type == 'GDT' & !control_tds$test) update_opt_st()
              if (control$Type == 'GDT' & control_tds$test) update_opt_st_test()
            } else if(control_tds$refine) {

              # --- Golden Search Refinement ---
              if (control$Type == 'GD'  ) {
                cat('\n one step Golden search ratio \n')
                if(k == tail(varying, 1)) converged = TRUE
                data$e0k <- data$e0 + BETA[, k] * X[, k]
                myformula_bk = as.formula(paste0('e0k~-1+', k))

                refined <- golden_search_bandwidth(
                  formula = myformula_bk, Ht = NULL, data = data, coords = coords,
                  fixed_vars = NULL, kernels = kernels, Model = "GWR", control = control,
                  lower.bound = HKmin[k], upper.bound = V[max(1, which(V == opt[k]) - 3)]
                )
                model_k <- refined$model
                opt[k] <- model_k@H
                betav <- model_k@Betav
                e0 <- residuals(model_k)
                isol <- is.na(betav)
                e0[isol] <- data$e0[isol]
                if (get_AIC) {
                  TSik[!isol, k] <- model_k@TS[!isol]
                  Sk <- model_k@Shat
                }
              }
              # (Note: GDT Refinement block commented out in original code, kept as such)
            } else update_opt_known()
          } else {
            update_opt_known()
          }

          # --- Update AIC trace and matrices ---
          if (control_tds$get_AIC) {
            Rkk[[k]] <- compute_Rk(Rk[[k]], Sk, St, foldsl)
            St = St + Rkk[[k]] - Rk[[k]]
            new_TS = diag(St)
            new_ts = sum(new_TS)

            if (verbose & control_tds$get_AIC)
              cat(' AICc : ', aicc_f(e0[idx_init], new_ts, n_time))
          }
          BETA[!isol, k] = betav[!isol]

          # --- Update residuals and RMSE ---
          fit = rowSums(BETA * X)
          data$e0 = Y - fit
          rmse = sqrt(mean(data$e0^2))
        }

        # ------------------------------------------------------------
        # 3.3 FIXED VARIABLES UPDATE
        # ------------------------------------------------------------
        for (k in fixed_vars) {
          last_rmse = rmse
          data$e0k <- data$e0 + BETA[, k] * X[, k]
          myformula_bk = as.formula(paste0('e0k~-1+', k))

          model_tds <- MGWRSAR(formula = myformula_bk, data = data,
                               coords = coords, fixed_vars = k, kernels = kernels,
                               H = NULL, Model = 'OLS', control = controlv)

          BETA[, k] = model_tds@Betac

          if (control_tds$get_AIC) {
            TSik[, k] <- unlist(res[mybest, 'TS'])
            Sk <- unlist(res[mybest, 'S'][[1]])
            Rkk[[k]] <- compute_Rk(Rk[[k]], Sk, St, foldsl)
            St = St + Rkk[[k]] - Rk[[k]]
            new_TS = diag(St)
            new_ts = sum(new_TS)
          }

          fit = rowSums(BETA * X)
          data$e0 = Y - fit
          rmse = sqrt(mean(data$e0^2))
        }

        # ------------------------------------------------------------
        # 3.4 UPDATE CRITERIA AND CONVERGENCE CHECKS
        # ------------------------------------------------------------
        if (control_tds$get_AIC)
          myAICc <- aicc_f(data$e0[idx_init], new_ts, n_time)

        # Update stability tracker
        if (control$Type == 'GDT') {
          stable[opt != last_opt | opt_t != last_opt_t] <- 0
          stable[opt == last_opt & opt_t == last_opt_t] <- stable[opt == last_opt & opt_t == last_opt_t] + 1
        } else {
          stable[opt != last_opt] <- 0
          stable[opt == last_opt] <- stable[opt == last_opt] + 1
        }

        # Compute deltas
        delta_rmse = (last_rmseG - sqrt(mean((Y - rowSums(BETA * X))^2))) / last_rmseG
        if (control_tds$get_AIC)
          delta_AICc = (last_AICc - myAICc) / last_AICc

        # ------------------------------------------------------------
        # 3.5 IF IMPROVEMENT: SAVE BEST STATE
        # ------------------------------------------------------------
        if(bestRMSE < sqrt(mean((Y - rowSums(BETA * X))^2))) patience=patience+1 else patience=0
        if (control_tds$get_AIC) {
          S = St
          for (k in varying) Rk[[k]] <- Rkk[[k]]
        }

        if (sqrt(mean((Y - rowSums(BETA * X))^2)) <= bestRMSE) {
          bestBETA = BETA
          bestRMSE = sqrt(mean((Y - rowSums(BETA * X))^2))
          H = opt
          if (control$Type == 'GDT') Ht = opt_t
          mybestG = i + 1
          if (control_tds$get_AIC) {
            mybestS = S
            mybestRk = Rk
          }
        }

        # ------------------------------------------------------------
        # 3.6 END-OF-ITERATION HOUSEKEEPING
        # ------------------------------------------------------------
        fit = rowSums(BETA * X)
        data$e0 = Y - fit
        rmse = sqrt(mean(data$e0^2))
        delta_rmse = (last_rmseG - rmse) / last_rmseG

        if (control_tds$get_AIC) {
          HTS[[i + 1]] <- new_TS
          HAICc <- c(HAICc, myAICc)
        }

        HBETA[[i + 1]] <- BETA

        if (!is.null(TRUEBETA)) {
          for (k in 1:K)
            HRMSE[i + 1, k] = sqrt(mean((TRUEBETA[, k] - BETA[, k])^2))
          HRMSE[i + 1, K + 1] <- mean(HRMSE[i + 1, 1:K])
          HRMSE[i + 1, K + 2] <- opt[k]
          if (control_tds$get_AIC) HRMSE[i + 1, K + 3] <- myAICc
          HRMSE[i + 1, K + 4] <- sqrt(mean(data$e0^2))
        }

        if (verbose & control_tds$get_AIC)
          cat('\n delta_AICc: ', delta_AICc, ' AICc : ', myAICc, ' stable ', stable, ' delta_rmse ', delta_rmse, ' rmse ', rmse)
        if (verbose & !control_tds$get_AIC)
          cat('\n stable ', stable, ' delta_rmse ', delta_rmse, ' rmse ', rmse)
        if (verbose) cat('\n\n H =', opt)
        if (verbose & control$Type == 'GDT') cat('\n Ht =', opt_t)

        # --- sign_history to detect oscillations ---
        current_sign <- sign(delta_rmse)
        if (!is.na(current_sign) && current_sign != 0) {
          sign_history <- c(sign_history, current_sign)
          delta_history <- c(delta_history, delta_rmse)
        }

        if (control$Type == "GDT")
          bandwidth_history[[length(bandwidth_history) + 1]] <- list(h = opt, ht = opt_t)
        else
          bandwidth_history[[length(bandwidth_history) + 1]] <- list(h = opt)

        # Limit history size
        if (length(sign_history) > 10) {
          sign_history <- tail(sign_history, 10)
          delta_history <- tail(delta_history, 10)
        }
        if (length(bandwidth_history) > 10)
          bandwidth_history <- tail(bandwidth_history, 10)

        # ============================================================
        # [1] Detect strict Yoyo  (RMSE + - + - +)
        # ============================================================
        if (length(sign_history) >= 4) {
          diffs <- diff(sign_history)
          alt_seq <- all(abs(diffs) == 2) && all(diff(diffs) != 0)
          alt_count <- if (alt_seq) length(diffs[diffs != 0]) else 0
        } else {
          alt_count <- 0
        }

        mean_last2 <- if (length(delta_history) >= 2) mean(abs(tail(delta_history, 2))) else abs(delta_rmse)

        # ============================================================
        # [2]  Detect ping-pong between bandwidths (repeated values)
        # ============================================================
        bw_pingpong <- 0
        if (length(bandwidth_history) >= 4) {
          for (j in seq_len(length(bandwidth_history) - 2)) {
            bw_prev <- bandwidth_history[[j + 1]]$h
            bw1 <- bandwidth_history[[j]]$h
            bw2 <- bandwidth_history[[j + 2]]$h

            has_change <- !isTRUE(all.equal(bw_prev, bw1, tolerance = 0))
            same_h <- isTRUE(all.equal(bw1, bw2, tolerance = 0))

            same_ht <- TRUE
            if (control$Type == "GDT") {
              bw_prev_t <- bandwidth_history[[j + 1]]$ht
              bw1_t <- bandwidth_history[[j]]$ht
              bw2_t <- bandwidth_history[[j + 2]]$ht
              has_change <- has_change || !isTRUE(all.equal(bw_prev_t, bw1_t, tolerance = 0))
              same_ht <- isTRUE(all.equal(bw1_t, bw2_t, tolerance = 0))
            }

            if (has_change && same_h && same_ht)
              bw_pingpong <- bw_pingpong + 1
          }
        }

        # ============================================================
        # [3] 3 times of perfect ping-pong -> INSERT INTERMEDIATE BW
        # ============================================================
        if (length(bandwidth_history) >= 3) {
          bw_tm2 <- bandwidth_history[[length(bandwidth_history) - 2L]]
          bw_tm1 <- bandwidth_history[[length(bandwidth_history) - 1L]]
          bw_t   <- bandwidth_history[[length(bandwidth_history)]]

          pp_spatial <- detect_perfect_pingpong(bw_tm2$h, bw_tm1$h, bw_t$h)
          pp_temporal <- TRUE
          if (control$Type == "GDT")
            pp_temporal <- detect_perfect_pingpong(bw_tm2$ht, bw_tm1$ht, bw_t$ht)

          if (pp_spatial && pp_temporal)
            pp_streak <- pp_streak + 1L
          else
            pp_streak <- 0L

          if (pp_streak >= 3L) {
            cat("\n[WARNING] Perfect ping-pong detected 3 times in a row -- inserting intermediate bandwidths.\n")

            # ----- SPATIAL CORRECTION -----
            h_now  <- as.numeric(bw_t$h)
            h_prev <- as.numeric(bw_tm1$h)

            if (isTRUE(control$adaptive[1])) {
              need_mid_h <- abs(h_now - h_prev) > 2
              h_mid <- round((h_now + h_prev) / 2)
            } else {
              rel_diff_h <- abs(h_now - h_prev) / pmax(1e-12, (abs(h_now) + abs(h_prev)) / 2)
              need_mid_h <- rel_diff_h > 0.001
              h_mid <- (h_now + h_prev) / 2
            }

            for (kk in seq_along(h_now)) {
              if (isTRUE(need_mid_h[kk])) {
                V <- sort(unique(c(V, h_mid[kk])))
                opt[kk] <- h_mid[kk]
              }
            }

            # ----- TEMPORAL CORRECTION -----
            if (control$Type == "GDT") {
              ht_now  <- as.numeric(bw_t$ht)
              ht_prev <- as.numeric(bw_tm1$ht)

              if (isTRUE(control$adaptive[2])) {
                need_mid_ht <- abs(ht_now - ht_prev) > 2
                ht_mid <- round((ht_now + ht_prev) / 2)
              } else {
                rel_diff_ht <- abs(ht_now - ht_prev) / pmax(1e-12, (abs(ht_now) + abs(ht_prev)) / 2)
                need_mid_ht <- rel_diff_ht > 0.001
                ht_mid <- (ht_now + ht_prev) / 2
              }

              for (kk in seq_along(ht_now)) {
                if (isTRUE(need_mid_ht[kk])) {
                  Vt <- sort(unique(c(Vt, ht_mid[kk])))
                  opt_t[kk] <- ht_mid[kk]
                }
              }
            }

            cat("# Intermediate bandwidths inserted. Resuming iterations.\n")
            pp_streak <- 0L
            sign_history <- integer(0)
            delta_history <- numeric(0)
            # Reset history to current intermediate state
            bandwidth_history <- list(list(h = opt, ht = if (control$Type == "GDT") opt_t else NULL))
          }
        }

        # ============================================================
        # [4] STOPPING CONDITIONS
        # ============================================================

        stop_flag <- FALSE
        reason <- NULL

        # (a) patience rule
        if (patience >4) {
          stop_flag <- TRUE
          reason <- sprintf("%.6f iteration withouts improvement",patience)
        }

        # (a) weak yoyo oscillation
        if (alt_count >= 3 && mean_last2 < tol) {
          stop_flag <- TRUE
          reason <- sprintf("weak yoyo oscillation (alt=%d, mean|RMSE|=%.6f < tol=%.6f)",
                            alt_count, mean_last2, tol)
        }

        # (b) repetition stricte des memes bandwidths (si la correction n'a pas suffi)
        if (bw_pingpong >= 3) {
          stop_flag <- TRUE
          reason <- sprintf("repeated bandwidth ping-pong (%d identical alternations)", bw_pingpong)
        }

        if (stop_flag) {
          cat("\n Stopping criterion reached:", reason, "\n")

          # Last resort average
          opt <- round((opt + last_opt) / 2)
          if (control$Type == "GDT")
            opt_t <- round((opt_t + last_opt_t) / 2)
          i = i + 1
          break
        }

        if (abs(delta_rmse) < tol) post_conv_count <- post_conv_count + 1 else post_conv_count <- 0
        i = i + 1
      }

      # ============================================================
      # [5] FINALIZATION & MODEL ASSEMBLY
      # ============================================================
      names(H) <- varying

      if (control_tds$get_AIC) {
        myAICc <- HAICc[mybestG - 1]
        TS = as.numeric(unlist(HTS[mybestG]))
        tS = sum(TS)
      }

      BETA <- bestBETA
      if (!is.null(TRUEBETA))
        HRMSE <- HRMSE[1:mybestG, ]
      HBETA <- HBETA[1:mybestG]

      fit = rowSums(BETA * X)
      data$e0 = Y - fit

      # 4.1 BUILD FINAL MODEL OBJECT
      modelGWR <- new('mgwrsar')
      modelGWR@mycall <- mycall
      modelGWR@data <- data
      modelGWR@coords <- as.matrix(coords)
      modelGWR@X <- modelGWR@XV <- X
      if (verbose) cat('\n Time = ', (proc.time() - start)[3], '\n')

      modelGWR@formula = formula
      modelGWR@Model = Model
      modelGWR@H = H
      modelGWR@fixed_vars <- unique(c(as.character(fixed_vars),
                                      names(modelGWR@H)[which(modelGWR@H == max_dist)]))
      modelGWR@Betav = BETA
      modelGWR@fit = rowSums(BETA * X)
      modelGWR@residuals = Y - modelGWR@fit
      modelGWR@RMSE = sqrt(mean(modelGWR@residuals^2))
      modelGWR@Type = control$Type
      modelGWR@kernels = kernels
      modelGWR@adaptive = control$adaptive
      modelGWR@TP = 1:n
      modelGWR@ncore = control_tds$ncore
      modelGWR@NN = control$NN

      if(control_tds$check_pairs) {
        if(!is.null(HKmin))  modelGWR@HKmin <- HKmin
        if(!is.null(HKMIN)) modelGWR@HKMIN <- HKMIN
      }
      if (!is.null(control$Z))
        modelGWR@Z = control$Z
      modelGWR@V = V

      if (control$Type == 'GDT') {
        modelGWR@Vt = Vt
        modelGWR@Ht = Ht
        modelGWR@alpha = control$alpha
      }

      modelGWR@Y = Y
      modelGWR@ctime <- (proc.time() - start)[3]
      modelGWR@HBETA <- HBETA

      if (control_tds$get_AIC) {
        if (!is.null(S)) modelGWR@Shat = mybestS
        modelGWR@TS = diag(mybestS)
        modelGWR@tS = sum(modelGWR@TS)
        modelGWR@edf <- n - tS
        modelGWR@AICc <- myAICc
        for (k in varying)
          modelGWR@R_k[[k]] <- mybestRk[[k]]
      }

      if (!is.null(TRUEBETA))
        modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[, 1]), ]

      if(Model == 'atds_mgwr') modelGWR@G = G
      control_tds$model_stage1 <- model_stage1 <- returned_model <- modelGWR
    })
  }

  myformula_b = as.formula(paste0('e0~-1+', paste0(varying, collapse = '+')))
  model_lm0 = model0
  model_stage1 = NULL
  control_tds$test = FALSE

  stage1_tds_mgwr()

  # ============================================================
  # 5. STAGE 2 - BOOSTING (OPTIONAL)
  # ============================================================
  stage2_atds_mgwr <- function(env = parent.frame()){
    with(env, {
      if(browser == -2) browser()
      if(verbose) cat('Start stage 2 : \n')
      nrounds = control_tds$nrounds
      if(!is.null(control_tds$model_stage1)) model_stage1 <- control_tds$model_stage1
      if(length(model_stage1@R_k) == 0) stop('Starting model_stage1 must be runned with get_AIC=TRUE')
      if(verbose) summary(model_stage1)
      control_tds$V <- model_stage1@V
      varying <- setdiff(namesX, fixed_vars)
      controlv <- control
      controlv$get_s = TRUE
      if(!is.null(control_tds$model_stage1)) {
        HBETA = model_stage1@HBETA
        HRMSE = model_stage1@HRMSE
        i <- nrow(HRMSE)
        HRMSE <- rbind(HRMSE, matrix(NA, ncol = ncol(HRMSE), nrow = nrounds + 1))
        BETA <- model_stage1@Betav
        AICc <- model_stage1@AICc
        varying <- varying[order(model_stage1@H, decreasing = T)]
        fixed_vars <- model_stage1@fixed_vars
        varyingT <- setdiff(varying, fixed_vars)
        data$e0 <- residuals(model_stage1)
        TS = model_stage1@TS
        H = model_stage1@H
        G <- model_stage1@G
        control$indexG = G$indexG
        control$dists = G$dists
      }
      if(verbose) cat('\n\n########################################################################\n STAGE 2 :  use multiple/adaptive bandwidth for each covariate \n using Top Down Scale Approach with backfitting ...\n########################################################################\n')
      HH <- list()

      nround = 1
      n <- nrow(model_stage1@Betav)
      if(control_tds$get_AIC) {
        Rkk <- Rk <- model_stage1@R_k
        S <- model_stage1@Shat
        TS <- model_stage1@TS
        tS = model_stage1@tS
        last_AICck <- last_AICc <- AICc + 1
      } else {
        AICc = 0;
        last_AICck <- last_AICc <- 1
        tS <- ts <- S <- NULL
      }
      rmse = sqrt(mean(data$e0^2)) + 10^6
      n_updated = 1


      while(nround <= nrounds & n_updated > 0) {
        last_AICck <- last_AICc <- AICc
        mybestbeta <- BETA
        mybestAICc <- AICc
        mybesttS <- tS
        mybestTS <- TS
        mybestShat = S
        mybestedf <- n - tS
        n_updated = 0
        for(k in namesX) {
          rmse <- sqrt(mean(data$e0^2))
          control_tds2 = control_tds
          control_tds2$H <- 0
          #control_tds2$H<-control_tds2$H[k]
          control_tds2$V <- control_tds2$V[-1]
          control_tds2$model_stage1 = NULL
          control_tds2$TRUEBETA = TRUEBETA
          if(!is.null(TRUEBETA)) {
            colnames(control_tds2$TRUEBETA) <- namesX
            control_tds2$TRUEBETA <- matrix(control_tds2$TRUEBETA[, k], ncol = 1)
          }
          data$e0k <- data$e0 + BETA[, k] * X[, k]
          myformula_bk = as.formula(paste0('e0k~-1+', k))
          if(k %in% fixed_vars) {
            model_tds <- MGWRSAR(formula = myformula_bk, data = data, coords = coords, fixed_vars = k, kernels = kernels, H = NULL, Model = 'OLS', control = controlv)
            HH[[k]] <- n
            model_tds@Betav <- as.matrix(rep(model_tds@Betac, n), ncol = 1)
          } else { model_tds <- atds_gwr(formula = myformula_bk, data = data, coords = coords, kernels = kernels, fixed_vars = NULL, control_tds = control_tds2, control = controlv)
          }
          betav <- model_tds@Betav
          e0 = residuals(model_tds)
          idxt = !is.na(betav)
          if(control_tds$get_AIC) {
            Sk <- model_tds@Shat
            Rkk[[k]] <- eigenMapMatMult(Sk, Rk[[k]]) + Sk - eigenMapMatMult(Sk, S)
            St = S + Rkk[[k]] - Rk[[k]]
            new_ts = sum(diag(St))
            model_tds@AICc <- n * log(sum(e0[idxt]^2) / n_time) + n_time * log(2 * pi) + n_time * (n_time + new_ts) / (n_time - 1 - new_ts)
          }
          if(nround < nrounds) fullupdate = FALSE else fullupdate = TRUE

          if( model_tds@AICc < last_AICck ) {
            if(verbose) cat(paste0(' ', k, ' updated; '))
            n_updated = n_updated + 1
            BETA[idxt, k] = betav[idxt]
            last_AICck <- AICc <- model_tds@AICc
            if(!(k %in% fixed_vars)) HH[[k]] <- model_tds@V[model_tds@V >= model_tds@H]


            fit = rowSums(BETA * X)
            data$e0 <- Y - fit
            if(control_tds$get_AIC) {
              S = St
              TS = diag(S)
              tS = sum(TS)
              Rk[[k]] <- Rkk[[k]]
            }
          } else {
            if(verbose)  cat(paste0(' ', k, ' not updated; '))
          }
        }
        if(nround == 1) {
          HM <- matrix(NA, nrow = length(V), ncol = length(varying))
          colnames(HM) <- varying

          for(k in varying) {
            if(!is.null(HH[[k]])) {
              minh <- min(unlist(HH[[k]]))
              vk <- V[V >= minh]
              HM[1:length(vk), k] <- vk
            }
          }
        }
        if(!(last_AICc - AICc) / abs(AICc) > tol & nrounds >= 10) {
          BETA <- mybestbeta
          AICc <- mybestAICc
          tS <- mybesttS
          TS <- mybestTS
          S <- mybestShat
        }
        HBETA[[length(HBETA) + 1]] <- BETA
        if(!is.null(TRUEBETA) ) {
          for(k in 1:K) HRMSE[i + 1, k] = sqrt(mean((TRUEBETA[, k] - BETA[, k])^2))
          HRMSE[i + 1, K + 1] <- mean(HRMSE[i + 1, 1:K])
          cat('\n BETA RMSE = ', mean(HRMSE[i + 1, 1:K]), '\n')
          HRMSE[i + 1, K + 2] <- min(unlist(HH))
          HRMSE[i + 1, K + 3] <- AICc
          HRMSE[i + 1, K + 4] <- sqrt(mean(data$e0^2))
        }
        nround = nround + 1
        i = i + 1
        if( !is.null(TRUEBETA) & verbose) cat('\n nround ', nround - 1, ' mean beta rmse ', HRMSE[i + 1, K + 1], '  rmse ', sqrt(mean(data$e0^2)), '\n')
        if( is.null(TRUEBETA)  & verbose ) cat('\n nround ', nround - 1,  ' rmse ', sqrt(mean(data$e0^2)), ' AICc ', AICc, '\n')
      }
      if(verbose) for(k in varying) {
        cat('\n', k, ' : ', unlist(HH[[k]] ), '\n')
      }
      ### RETURN MODEL
      modelGWR <- model_stage1
      modelGWR@mycall <- mycall
      modelGWR@data <- data
      modelGWR@coords <- as.matrix(coords)
      modelGWR@X <- modelGWR@XV <- X
      modelGWR@formula = formula
      if(verbose)  cat('\n Time = ', (proc.time() - start)[3], '\n')
      modelGWR@Model = Model
      modelGWR@fixed_vars <- as.character(fixed_vars)
      modelGWR@Betav = BETA
      modelGWR@fit = rowSums(BETA * X)
      modelGWR@residuals = Y - modelGWR@fit
      modelGWR@RMSE = sqrt(mean(modelGWR@residuals^2))


      modelGWR@HM = HM
      modelGWR@H = model_stage1@H
      if(model_stage1@Type == 'GDT') {
        modelGWR@Type = model_stage1@Type
        modelGWR@kernels = model_stage1@kernels
        modelGWR@Ht = model_stage1@Ht
        modelGWR@adaptive = model_stage1@adaptive
      } else {
        modelGWR@Type = control$Type
        modelGWR@kernels = kernels
        modelGWR@adaptive = control$adaptive
      }
      modelGWR@TP = 1:n
      modelGWR@ncore = control_tds$ncore
      modelGWR@TP = 1:n
      modelGWR@V = c(n, V)
      modelGWR@Y = Y
      modelGWR@ctime <-  (proc.time() - start)[3]
      modelGWR@HBETA <- HBETA
      if(!is.null(TRUEBETA)) {
        modelGWR@HRMSE <- HRMSE[!is.na(HRMSE[, 1]), ]
      }
      if(control_tds$get_AIC) {
        modelGWR@AICc <- AICc
        modelGWR@tS = tS
        modelGWR@TS = TS
        modelGWR@Shat = S
        modelGWR@edf <- n - tS
      }
      returned_model <- modelGWR
    })
  }

  # ============================================================
  # TRANSITION & STAGE 2 EXECUTION
  # ============================================================

  if (Model == 'atds_mgwr') {
    message("Running Stage 2 (ATDS-MGWR boosting)...")

    # 1. On passe le modèle du stage 1 dans les paramètres de contrôle pour le stage 2
    # (car stage2_atds_mgwr cherche 'model_stage1' dans control_tds)
    control_tds$model_stage1 <- returned_model

    # 2. On appelle la bonne fonction (mgwr pas gwr) sans arguments (utilise parent.frame)
    stage2_atds_mgwr()
  }

  gc()
  return(returned_model)
}
