#' Optimization of 2D Bandwidths (Spatial and Temporal) using Golden Section Search
#'
#' @description
#' This function optimizes spatial and temporal bandwidths simultaneously using
#' a 2D Golden Section Search approach. It is typically used internally by
#' \code{search_bandwidths} when \code{refine = TRUE} for
#' spatio-temporal models (Type 'GDT').
#'
#' @usage golden_search_2d_bandwidth(formula, data, coords, fixed_vars, kernels, Model,
#'                                   control, lower.bound.space, upper.bound.space,
#'                                   lower.bound.time, upper.bound.time,
#'                                   tolerance_s = 1e-06, tolerance_t = 1e-06,
#'                                   max_iter = 10)
#'
#' @param formula A formula object.
#' @param data A data frame containing the variables.
#' @param coords A matrix of coordinates (spatial indices).
#' @param fixed_vars Vector of names of variables with fixed coefficients.
#' @param kernels Vector of kernel types (e.g., c('gauss', 'gauss')).
#' @param Model Character string specifying the model type (e.g., 'GWR').
#' @param control List of control parameters.
#' @param lower.bound.space Numeric. Lower bound for the spatial bandwidth search.
#' @param upper.bound.space Numeric. Upper bound for the spatial bandwidth search.
#' @param lower.bound.time Numeric. Lower bound for the temporal bandwidth search.
#' @param upper.bound.time Numeric. Upper bound for the temporal bandwidth search.
#' @param tolerance_s Numeric. Convergence tolerance for the spatial dimension. Default is 1e-6.
#' @param tolerance_t Numeric. Convergence tolerance for the temporal dimension. Default is 1e-6.
#' @param max_iter Integer. Maximum number of iterations for the optimization loop. Default is 10.
#'
#' @return A list containing:
#' \item{minimum}{Vector of optimized bandwidths c(h_spatial, h_temporal).}
#' \item{objective}{Final AICc value.}
#' \item{model}{The fitted mgwrsar model object.}
#'
#' @seealso \code{\link{golden_search_bandwidth}}, \code{\link{search_bandwidths}}
#' @export
golden_search_2d_bandwidth <- function(
    formula, data, coords, fixed_vars, kernels,
    Model, control,
    lower.bound.space, upper.bound.space,
    lower.bound.time,  upper.bound.time,
    tolerance_s = 1e-6, tolerance_t = 1e-6,
    max_iter = 10
) {

  mf <- model.frame(formula, data)
  data <- data[,names(mf)]

  # -----------------------------------------------------------
  # Rounding helpers: impose tolerance-based bandwidth precision
  # -----------------------------------------------------------
  round_space <- function(x) round(x / tolerance_s) * tolerance_s
  round_time  <- function(x) round(x / tolerance_t) * tolerance_t

  golden_ratio <- 2/(sqrt(5)+1)

  # Slight jitter if duplicated coords exist
  set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
  if(control$Type %in% c('GD','GDT')) if(sum(duplicated(coords))>0) {
    set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
    coords[,1] <- coords[,1] + seq_along(coords[,1]) * 1e-12
    coords[,2] <- coords[,2] + seq_along(coords[,2]) * 1e-12
  }
  if(control$Type %in% c('GDT','T')) if(sum(duplicated(control$Z))>0) {
    set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
    control$Z<- control$Z + seq_along(control$Z) * 1e-12
  }

  # Adaptive spec
  if (length(control$adaptive) != 2)
    stop("control$adaptive must be c(adaptive_space, adaptive_time).")

  adaptive_space <- control$adaptive[1]
  adaptive_time  <- control$adaptive[2]

  # -----------------------------------------------------------
  # Lazy storage for corner cases (GD / T)
  # -----------------------------------------------------------
  corner_controls <- list(
    controlvd = NULL,
    controlvt = NULL,
    dd_GD     = NULL,
    dd_T      = NULL
  )

  # -----------------------------------------------------------
  # AICc with post-estimation penalty (variance + collinearity)
  # -----------------------------------------------------------
  AICc_with_penalty <- function(model, control) {

    e   <- model@residuals
    RSS <- sum(e^2)
    n   <- length(e)
    trS <- model@tS

    # Raw AICc formula
    AICc_raw <- n*log(RSS/n) + n*(n+trS)/(n-trS-2)

    # ---------------------------------------------------
    # Penalties (if rules are provided in control)
    # ---------------------------------------------------
    pen <- 0

    # --- [A] Variance local penalty ---
    if (!is.null(control$pen_check_var) && control$pen_check_var) {
      ok_var <- check_local_variance_post(
        model = model,
        X     = model@X,          # MGWRSAR stores X inside
        var_ratio_min = control$min_var_ratio  %||% 0.05,
        var_q         = control$min_var_q      %||% 0.10,
        n_loc_sample  = control$n_loc_sample   %||% 30L
      )
      if (!ok_var)
        pen <- pen + (control$lambda_var %||% 1e6)
    }

    # --- [B] Local collinearity penalty ---
    if (!is.null(control$pen_check_kappa) && control$pen_check_kappa) {
      ok_kappa <- check_local_condition_post(
        model = model,
        X     = model@X,
        kappa_max = control$max_kappa %||% 1e8,
        n_loc_sample = control$n_loc_sample %||% 30L
      )
      if (!ok_kappa)
        pen <- pen + (control$lambda_kappa %||% 1e6)
    }

    # --- [C] Trace(S) penalty ---
    if (!is.null(control$max_trS)) {
      if (trS > control$max_trS) {
        pen <- pen + (control$lambda_trS %||% 1e6)*(trS - control$max_trS)^2
      }
    }

    return(AICc_raw + pen)
  }

  # -----------------------------------------------------------
  # FAST post-MGWRSAR VARIANCE CHECK (multivariate)
  # -----------------------------------------------------------
  check_local_variance_post <- function(model, X,
                                        var_ratio_min = 0.05,
                                        var_q = 0.10,
                                        n_loc_sample = 30L) {

    Wi <- model@Wi    # list of weight vectors per location
    if (is.null(Wi)) return(TRUE)

    var_glob <- apply(X, 2, var)
    if (any(!is.finite(var_glob) | var_glob <= .Machine$double.eps))
      return(FALSE)

    n <- length(Wi)
    idx <- if (n <= n_loc_sample) seq_len(n) else sample.int(n, n_loc_sample)

    for (i in idx) {
      w <- Wi[[i]]
      if (sum(w) < 1e-12) next

      wx  <- drop(t(w) %*% X) / sum(w)
      wx2 <- drop(t(w) %*% (X^2)) / sum(w)
      ratio <- (wx2 - wx^2) / var_glob

      qv <- quantile(ratio, var_q, na.rm = TRUE)
      if (qv < var_ratio_min) return(FALSE)
    }

    TRUE
  }

  # -----------------------------------------------------------
  # FAST post-MGWRSAR COLLINEARITY CHECK via X'WX
  # -----------------------------------------------------------
  check_local_condition_post <- function(model, X,
                                         kappa_max = 1e8,
                                         n_loc_sample = 30L) {

    Wi <- model@Wi
    if (is.null(Wi)) return(TRUE)

    n  <- length(Wi)
    p  <- ncol(X)
    idx <- if (n <= n_loc_sample) seq_len(n) else sample.int(n, n_loc_sample)

    for (i in idx) {
      w <- Wi[[i]]
      if (sum(w) < 1e-12) next
      Xw <- sqrt(w) * X
      XtWX <- crossprod(Xw)

      k_loc <- tryCatch(kappa(XtWX, exact = FALSE), error = function(e) Inf)
      if (!is.finite(k_loc) || k_loc > kappa_max)
        return(FALSE)
    }

    TRUE
  }

  # ==========================================================
  # GOLDEN SEARCH IN 1D (inner loop)
  # ==========================================================
  golden_section_1d <- function(
    fixed_h, lower, upper, adaptive, is_space,
    formula, data, coords, fixed_vars, kernels, Model, control,
    tolerance, golden_ratio,
    lower.bound.space, upper.bound.space,
    lower.bound.time,  upper.bound.time
  ) {

    # Initial golden points
    x1 <- upper - golden_ratio*(upper - lower)
    x2 <- lower + golden_ratio*(upper - lower)

    if (adaptive) {
      x1 <- floor(x1); x2 <- ceiling(x2)
    } else {
      if (is_space) { x1 <- round_space(x1); x2 <- round_space(x2) }
      else          { x1 <- round_time(x1);  x2 <- round_time(x2) }
    }

    # ---------------------------------------
    # eval_cv : local AICc + penalties
    # ---------------------------------------
    eval_cv <- function(h1, fixed_h) {

      # Rebuild hs, ht
      if (is_space) { hs <- h1; ht <- fixed_h }
      else          { hs <- fixed_h; ht <- h1 }

      # Determine which control to use (full / GD / T / OLS)
      # FULL (GDT)
      if (hs < upper.bound.space && ht < upper.bound.time) {
        ctrl_use    <- control
        kernels_use <- kernels
        Model_use   <- "GWR"
      }
      # OLS corner
      else if (hs >= upper.bound.space && ht >= upper.bound.time) {
        ctrl_use    <- control
        kernels_use <- kernels
        Model_use   <- "OLS"
      }
      # GD corner
      else if (hs < upper.bound.space && ht >= upper.bound.time) {
        if (is.null(corner_controls$controlvd)) {
          ctrl <- modifyList(control, list(
            dists   = NULL, indexG = NULL,
            Type    = "GD",
            adaptive = control$adaptive[1]
          ))
          dd <- prep_d(coords  = coords,
                       NN      = control$NN,
                       TP      = control$TP,
                       kernels = kernels[1],
                       Type    = "GD")
          ctrl$dists  <- dd$dists
          ctrl$indexG <- dd$indexG
          corner_controls$controlvd <- ctrl
          corner_controls$dd_GD     <- dd
        }
        ctrl_use    <- corner_controls$controlvd
        kernels_use <- kernels
        Model_use   <- "GWR"
      }
      # T corner
      else {
        if (is.null(corner_controls$controlvt)) {
          ctrl <- modifyList(control, list(
            dists   = NULL, indexG = NULL,
            Type    = "T",
            adaptive = FALSE
          ))
          dd <- prep_d(coords  = as.matrix(control$Z, ncol=1),
                       NN      = control$NN,
                       TP      = seq_len(nrow(coords)),
                       kernels = kernels[2],
                       Type    = "T")
          ctrl$dists  <- dd$dists
          ctrl$indexG <- dd$indexG
          corner_controls$controlvt <- ctrl
          corner_controls$dd_T      <- dd
        }
        ctrl_use    <- corner_controls$controlvt
        kernels_use <- kernels
        Model_use   <- "GWR"
      }

      # ----------------------------------------------------
      # FULL MGWRSAR estimation for this (hs, ht)
      # ----------------------------------------------------
      model_loc <- tryCatch(
        MGWRSAR(
          formula, data, coords, fixed_vars,
          kernels_use,
          H = c(hs, ht),
          Model = Model_use, control = ctrl_use
        ),
        error = function(e) NULL
      )
      if (is.null(model_loc)) return(Inf)

      # ----------------------------------------------------
      # AICc with penalties (variance, kappa, trS)
      # ----------------------------------------------------
      return(AICc_with_penalty(model_loc, ctrl_use))
    }
    # Initial evaluations
    f1 <- eval_cv(x1, fixed_h)
    f2 <- eval_cv(x2, fixed_h)

    # Golden loop
    while ((abs(upper - lower) > tolerance) &&
           (!adaptive || abs(x2 - x1) > 1)) {

      if (f2 > f1) {
        upper <- x2
        x2 <- x1
        f2 <- f1
        x1 <- upper - golden_ratio*(upper - lower)
        if (adaptive) x1 <- floor(x1)
        else          x1 <- if (is_space) round_space(x1) else round_time(x1)
        f1 <- eval_cv(x1, fixed_h)
      } else {
        lower <- x1
        x1 <- x2
        f1 <- f2
        x2 <- lower + golden_ratio*(upper - lower)
        if (adaptive) x2 <- ceiling(x2)
        else          x2 <- if (is_space) round_space(x2) else round_time(x2)
        f2 <- eval_cv(x2, fixed_h)
      }
    }

    # Final selection
    if (adaptive)
      h_opt <- if (f1 < f2) x1 else x2
    else {
      h_opt <- (lower + upper)/2
      h_opt <- if (is_space) round_space(h_opt) else round_time(h_opt)
    }

    list(h_opt = h_opt, lower = lower, upper = upper)
  } # end golden_section_1d



  # ==========================================================
  # PRECOMPUTED FULL ST DISTANCES IF NEEDED
  # ==========================================================
  if (is.null(control$dists)) {

    if (is.null(control$NN)) control$NN <- nrow(data)
    if (is.null(control$TP)) control$TP <- seq_len(nrow(data))

    S <- if (length(kernels) > 1)
      as.matrix(cbind(coords, control$Z))
    else
      as.matrix(coords)

    stage1 <- prep_d(
      coords = S,
      NN = control$NN, TP = control$TP,
      kernels = kernels,
      Type = control$Type
    )

    control$dists  <- stage1$dists
    control$indexG <- stage1$indexG
  }


  # ==========================================================
  # OUTER 2D COORDINATEWISE GOLDEN SEARCH
  # ==========================================================
  h_space <- (lower.bound.space + upper.bound.space)/2
  h_time  <- (lower.bound.time  + upper.bound.time)/2

  for (iteration in 1:max_iter) {

    # ---- 1. Spatial update ----
    opt_space <- golden_section_1d(
      fixed_h = h_time,
      lower = lower.bound.space,
      upper = upper.bound.space,
      adaptive = adaptive_space,
      is_space = TRUE,
      formula, data, coords, fixed_vars, kernels,
      Model, control,
      tolerance_s, golden_ratio,
      lower.bound.space, upper.bound.space,
      lower.bound.time,  upper.bound.time
    )

    h_space_new <- opt_space$h_opt
    lower.bound.space <- (lower.bound.space + opt_space$lower)/2
    upper.bound.space <- (upper.bound.space + opt_space$upper)/2

    # ---- 2. Temporal update ----
    opt_time <- golden_section_1d(
      fixed_h = h_space_new,
      lower = lower.bound.time,
      upper = upper.bound.time,
      adaptive = adaptive_time,
      is_space = FALSE,
      formula, data, coords, fixed_vars, kernels,
      Model, control,
      tolerance_t, golden_ratio,
      lower.bound.space, upper.bound.space,
      lower.bound.time,  upper.bound.time
    )

    h_time_new <- opt_time$h_opt
    lower.bound.time <- (lower.bound.time + opt_time$lower)/2
    upper.bound.time <- (upper.bound.time + opt_time$upper)/2

    if (isTRUE(control$verbose)) {
      cat(sprintf(
        "Iteration %d: h_space = %.4f, h_time = %.4f\n",
        iteration, h_space_new, h_time_new
      ))
    }

    # Convergence check
    conv_space <- if (adaptive_space)
      abs(h_space_new - h_space) < 1
    else
      abs(h_space_new - h_space) < tolerance_s

    conv_time <- if (adaptive_time)
      abs(h_time_new - h_time) < 1
    else
      abs(h_time_new - h_time) < tolerance_t

    h_space <- h_space_new
    h_time  <- h_time_new

    if (conv_space && conv_time) break
  }


  # ==========================================================
  # FINAL MODEL
  # ==========================================================
  if (!is.finite(h_space) || !is.finite(h_time))
    stop("Golden search failed: non-finite bandwidths.")

  model <- MGWRSAR(
    formula, data, coords, fixed_vars, kernels,
    H = c(h_space, h_time),
    Model = Model, control = control
  )

  model@kernels <- kernels
  model@Z       <- control$Z
  model@mycall  <- match.call()

  list(
    minimum = c(h_space, h_time),
    model   = model
  )
}

golden_search_2d_bandwidth_old0<- function(
    formula, data, coords, fixed_vars, kernels,
    Model, control,
    lower.bound.space, upper.bound.space,
    lower.bound.time,  upper.bound.time,
    tolerance_s = 1e-6, tolerance_t = 1e-6,
    max_iter = 10
) {

  # --------------------------------------
  # Rounding helpers based on tolerances
  # --------------------------------------
  round_space <- function(x) {
    round(x / tolerance_s) * tolerance_s
  }
  round_time <- function(x) {
    round(x / tolerance_t) * tolerance_t
  }

  Z <- control$Z
  golden_ratio <- 2 / (sqrt(5) + 1)

  # Avoid duplicates in coords
  if (sum(duplicated(coords)) > 0)
    coords <- jitter(coords, amount = 0.0000001)

  # Proper adaptive vector
  if (length(control$adaptive) != 2)
    stop("control$adaptive must be c(adaptive_space, adaptive_time).")

  adaptive_space <- control$adaptive[1]
  adaptive_time  <- control$adaptive[2]

  # Lazy-memory containers for corner controls
  corner_controls <- list(
    controlvd = NULL,  # GD
    controlvt = NULL,  # T
    dd_GD     = NULL,
    dd_T      = NULL
  )

  # ==========================================================
  # GOLDEN SECTION SEARCH (1D)
  # ==========================================================
  golden_section_1d <- function(
    fixed_h, lower, upper, adaptive, is_space,
    formula, data, coords, fixed_vars, kernels, Model, control,
    tolerance, golden_ratio,
    lower.bound.space, upper.bound.space,
    lower.bound.time,  upper.bound.time
  ) {

    # initial golden points
    x1 <- upper - golden_ratio * (upper - lower)
    x2 <- lower + golden_ratio * (upper - lower)

    # adaptive mode: discrete NN / TP
    if (adaptive) {
      x1 <- floor(x1); x2 <- ceiling(x2)
    } else {
      # rounding candidate bandwidths
      if (is_space) {
        x1 <- round_space(x1); x2 <- round_space(x2)
      } else {
        x1 <- round_time(x1); x2 <- round_time(x2)
      }
    }
    # ======================================================
    # INTERNAL CV evaluation with corner logic
    # ======================================================
    eval_cv <- function(h1, fixed_h) {

      # reconstruct (hs, ht)
      if (is_space) {
        hs <- h1
        ht <- fixed_h
      } else {
        hs <- fixed_h
        ht <- h1
      }

      # FULL GWR
      if (hs < upper.bound.space && ht < upper.bound.time) {
        return(AICc_CV(
          c(hs, ht), formula, data, coords,
          fixed_vars, kernels, "GWR", control
        ))
      }

      # OLS corner
      if (hs >= upper.bound.space && ht >= upper.bound.time) {
        return(AICc_CV(
          c(hs, ht), formula, data, coords,
          fixed_vars, kernels, "OLS", control
        ))
      }

      # -------- GD corner (spatial-only) --------
      if (hs < upper.bound.space && ht >= upper.bound.time) {

        if (is.null(corner_controls$controlvd)) {

          ctrl <- modifyList(control, list(
            dists   = NULL,
            indexG  = NULL,
            Type    = "GD",
            adaptive = control$adaptive[1]
          ))

          dd <- prep_d(
            coords  = coords,
            NN      = control$NN,
            TP      = control$TP,
            kernels = kernels[1],
            Type    = "GD"
          )

          ctrl$dists  <- dd$dists
          ctrl$indexG <- dd$indexG

          corner_controls$controlvd <- ctrl
          corner_controls$dd_GD <- dd
        }

        return(AICc_CV(
          c(hs, ht), formula, data, coords,
          fixed_vars, kernels, "GWR",
          corner_controls$controlvd
        ))
      }

      # -------- T corner (temporal-only) --------
      if (hs >= upper.bound.space && ht < upper.bound.time) {

        if (is.null(corner_controls$controlvt)) {

          ctrl <- modifyList(control, list(
            dists   = NULL,
            indexG  = NULL,
            Type    = "T",
            adaptive = FALSE
          ))

          dd <- prep_d(
            coords  = as.matrix(control$Z, ncol = 1),
            NN      = control$NN,
            TP      = 1:nrow(coords),
            kernels = kernels[2],
            Type    = "T"
          )

          ctrl$dists  <- dd$dists
          ctrl$indexG <- dd$indexG

          corner_controls$controlvt <- ctrl
          corner_controls$dd_T <- dd
        }

        return(AICc_CV(
          c(hs, ht), formula, data, coords,
          fixed_vars, kernels, "GWR",
          corner_controls$controlvt
        ))
      }

      return(Inf)
    }

    # initial evaluations
    f1 <- eval_cv(x1, fixed_h)
    f2 <- eval_cv(x2, fixed_h)

    # ============================
    # Golden Iteration Loop
    # ============================
    while ((abs(upper - lower) > tolerance) &&
           (!adaptive || abs(x2 - x1) > 1)) {

      if (f2 > f1) {
        upper <- x2
        x2 <- x1
        f2 <- f1
        x1 <- upper - golden_ratio * (upper - lower)

        if (adaptive) { x1 <- floor(x1)
        } else {
          x1 <- if (is_space) round_space(x1) else round_time(x1)
        }

        f1 <- eval_cv(x1, fixed_h)

      } else {

        lower <- x1
        x1 <- x2
        f1 <- f2
        x2 <- lower + golden_ratio * (upper - lower)

        if (adaptive) { x2 <- ceiling(x2)
        } else {
          x2 <- if (is_space) round_space(x2) else round_time(x2)
        }

        f2 <- eval_cv(x2, fixed_h)
      }
    }

    # final h_opt
    if (adaptive) {
      h_opt <- if (f1 < f2) x1 else x2
    } else {
      h_opt <- (lower + upper) / 2
      h_opt <- if (is_space) round_space(h_opt) else round_time(h_opt)
    }

    list(h_opt = h_opt, lower = lower, upper = upper)
  }
  # ======================================================
  # Full spatio-temporal distances if not provided
  # ======================================================
  if (is.null(control$dists)) {

    if (is.null(control$NN)) control$NN <- nrow(data)
    if (is.null(control$TP)) control$TP <- 1:nrow(data)

    S <- if (length(kernels) > 1)
      as.matrix(cbind(coords, control$Z))
    else
      as.matrix(coords)

    stage1 <- prep_d(
      coords = S, NN = control$NN, TP = control$TP,
      kernels = kernels, Type = control$Type
    )
    control$indexG <- stage1$indexG
    control$dists  <- stage1$dists
  }

  # ======================================================
  # 2D ALTERNATING GOLDEN SEARCH
  # ======================================================
  h_space <- (lower.bound.space + upper.bound.space) / 2
  h_time  <- (lower.bound.time  + upper.bound.time)  / 2

  for (iteration in 1:max_iter) {

    # --- optimize spatial bandwidth ---
    opt_space <- golden_section_1d(
      fixed_h = h_time,
      lower = lower.bound.space, upper = upper.bound.space,
      adaptive = adaptive_space, is_space = TRUE,
      formula, data, coords, fixed_vars, kernels, Model, control,
      tolerance_s, golden_ratio,
      lower.bound.space, upper.bound.space,
      lower.bound.time,  upper.bound.time
    )

    h_space_new <- opt_space$h_opt
    lower.bound.space <- (lower.bound.space + opt_space$lower) / 2
    upper.bound.space <- (upper.bound.space + opt_space$upper) / 2

    # --- optimize temporal bandwidth ---
    opt_time <- golden_section_1d(
      fixed_h = h_space_new,
      lower = lower.bound.time, upper = upper.bound.time,
      adaptive = adaptive_time, is_space = FALSE,
      formula, data, coords, fixed_vars, kernels, Model, control,
      tolerance_t, golden_ratio,
      lower.bound.space, upper.bound.space,
      lower.bound.time,  upper.bound.time
    )

    h_time_new <- opt_time$h_opt
    lower.bound.time <- (lower.bound.time + opt_time$lower) / 2
    upper.bound.time <- (upper.bound.time + opt_time$upper) / 2

    if (control$verbose) {
      cat(sprintf("Iteration %d: h_space=%.4f, h_time=%.4f\n",
                  iteration, h_space_new, h_time_new))
    }

    conv_space <- if (adaptive_space)
      abs(h_space_new - h_space) < 1
    else
      abs(h_space_new - h_space) < tolerance_s

    conv_time <- if (adaptive_time)
      abs(h_time_new - h_time) < 1
    else
      abs(h_time_new - h_time) < tolerance_t

    h_space <- h_space_new
    h_time  <- h_time_new

    if (conv_space && conv_time) break
  }

  # final checks and model
  if (!is.finite(h_space) || !is.finite(h_time))
    stop("Optimization failed: non-finite bandwidths.")

  model <- MGWRSAR(
    formula, data, coords, fixed_vars, kernels,
    H = c(h_space, h_time), Model = Model, control = control
  )

  model@kernels <- kernels
  model@Z       <- Z
  model@mycall  <- match.call()

  list(
    minimum = c(h_space, h_time),
    model   = model
  )
}


golden_search_2d_bandwidth_old <- function(formula, data, coords, fixed_vars, kernels, Model, control,
                                       lower.bound.space, upper.bound.space,
                                       lower.bound.time, upper.bound.time,
                                       tolerance_s = 1e-6, tolerance_t = 1e-6, max_iter = 10) {
  Z<-control$Z
  golden_section_1d <- function(fixed_h, lower, upper, adaptive, is_space,
                                formula, data, coords, fixed_vars, kernels, Model, control,
                                tolerance, golden_ratio) {

    x1 <- upper - golden_ratio * (upper - lower)
    x2 <- lower + golden_ratio * (upper - lower)

    if (adaptive) {
      x1 <- floor(x1)
      x2 <- ceiling(x2)
    }

    eval_cv <- function(h1, h2,upper.bound.space,upper.bound.time) {
      tryCatch(
        ### a faire introduire les cas T et GD
        if (v < upper.bound.space & vt < upper.bound.time) {
          AICc_CV(
            if (is_space) c(h1, fixed_h) else c(fixed_h, h1),
            formula, data, coords, fixed_vars = NULL, kernels, Model = "GWR", control
          )
        } else if (v == upper.bound.space & vt == upper.bound.time) {
          AICc_CV(
            if (is_space) c(h1, fixed_h) else c(fixed_h, h1),
            formula, data, coords, fixed_vars = NULL, kernels, Model = "OLS", control
          )
        } else if (v < upper.bound.space & vt == upper.bound.time) {
          if (!exists("controlvd", envir =  parent.frame(), inherits = FALSE)) {
            controlvd <- modifyList(controlv, list(
              dists = NULL, indexG = NULL, Type = "GD",
              adaptive = controlv$adaptive[1]
            ))
          }
          AICc_CV(
            if (is_space) c(h1, fixed_h) else c(fixed_h, h1),
            formula, data, coords, fixed_vars = NULL, kernels, Model = "GWR", controlvd
          )
        } else if (v == upper.bound.space & vt < upper.bound.time) {
          if (!exists("controlvt",  envir =  parent.frame(),inherits = FALSE)) {
            controlvt <- modifyList(controlv, list(
              dists = NULL, indexG = NULL, Type = "T", adaptive = FALSE
            ))
          }
          AICc_CV(
            if (is_space) c(h1, fixed_h) else c(fixed_h, h1),
            formula, data, coords, fixed_vars = NULL, kernels, Model = "GWR", controlvt
          )
        }
        # AICc_CV(
        #   if (is_space) c(h1, fixed_h) else c(fixed_h, h1),
        #   formula, data, coords, fixed_vars, kernels, Model, control
        # )
        ,
        error = function(e) Inf
      )
    }

    f1 <- eval_cv(x1, fixed_h)
    f2 <- eval_cv(x2, fixed_h)

    while ((abs(upper - lower) > tolerance) &&
           (!adaptive || abs(x2 - x1) > 1)) {

      if (f2 > f1) {
        upper <- x2
        x2 <- x1
        f2 <- f1
        x1 <- upper - golden_ratio * (upper - lower)
        if (adaptive) x1 <- floor(x1)
        f1 <- eval_cv(x1, fixed_h)
      } else {
        lower <- x1
        x1 <- x2
        f1 <- f2
        x2 <- lower + golden_ratio * (upper - lower)
        if (adaptive) x2 <- ceiling(x2)
        f2 <- eval_cv(x2, fixed_h)
      }
    }

    h_opt <- if (adaptive) if (f1 < f2) x1 else x2 else (lower + upper) / 2
    list(h_opt = h_opt, lower = lower, upper = upper)
  }

  ptm <- proc.time()
  set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")

  if (sum(duplicated(coords)) > 0) {
    coords <- jitter(coords, amount = 0.0000001)
  }

  if (length(control$adaptive) != 2) {
    stop("control$adaptive must have length 2: [spatial, temporal]")
  }
  if (is.null(control$Z)) {
    stop("control$Z must be provided for temporal coordinates.")
  }

  adaptive_space <- control$adaptive[1]
  adaptive_time  <- control$adaptive[2]
  golden_ratio   <- 2 / (sqrt(5) + 1)

  h_space <- (lower.bound.space + upper.bound.space) / 2
  h_time  <- (lower.bound.time + upper.bound.time) / 2


  ## distance computation
  if(is.null(control$dists)){
    if(is.null(control$NN)) control$NN=nrow(data)
    if(is.null(control$TP)) control$TP=1:nrow(data)
    if (length(kernels) > 1) S = as.matrix(cbind(coords, control$Z)) else S = as.matrix(coords)
    stage1=prep_d(coords=S,NN=control$NN,TP=control$TP,kernels=kernels,Type=control$Type)
    control$indexG=stage1$indexG
    control$dists=stage1$dists
  }

  for (iteration in 1:max_iter) {
    #### 1. Optimize Spatial Bandwidth (fixing h_time) ####
    opt_space <- golden_section_1d(
      fixed_h = h_time,
      lower = lower.bound.space, upper = upper.bound.space,
      adaptive = adaptive_space, is_space = TRUE,
      formula = formula, data = data, coords = coords,
      fixed_vars = fixed_vars, kernels = kernels, Model = Model, control = control,
      tolerance = tolerance_s, golden_ratio = golden_ratio
    )

    h_space_new <- opt_space$h_opt
    lower.bound.space <- (lower.bound.space + opt_space$lower)/2
    upper.bound.space <-  (upper.bound.space+opt_space$upper)/2

    #### 2. Optimize Temporal Bandwidth (fixing h_space) ####
    opt_time <- golden_section_1d(
      fixed_h = h_space_new,
      lower = lower.bound.time, upper = upper.bound.time,
      adaptive = adaptive_time, is_space = FALSE,
      formula = formula, data = data, coords = coords,
      fixed_vars = fixed_vars, kernels = kernels, Model = Model, control = control,
      tolerance = tolerance_t, golden_ratio = golden_ratio
    )

    h_time_new <- opt_time$h_opt
    lower.bound.time <- (lower.bound.time+opt_time$lower)/2
    upper.bound.time <- (upper.bound.time +opt_time$upper)/2

    if (control$verbose) {
      cat(sprintf("Iteration %d: h_space = %.4f, h_time = %.4f\n", iteration, h_space_new, h_time_new))
    }

    # Check convergence
    conv_space <- if (adaptive_space) abs(h_space_new - h_space) < 1 else abs(h_space_new - h_space) < tolerance_s
    conv_time  <- if (adaptive_time) abs(h_time_new - h_time) < 1 else abs(h_time_new - h_time) < tolerance_t

    h_space <- h_space_new
    h_time  <- h_time_new

    if (conv_space && conv_time) break
  }

  # Final model
  if (!is.finite(h_space) || !is.finite(h_time)) {
    stop("Optimization failed: bandwidths are not finite.")
  }

  model <- MGWRSAR(formula, data, coords, fixed_vars, kernels, H = c(h_space, h_time), Model = Model, control = control)
  model@kernels<-kernels
  model@Z<-Z
  model@mycall <- match.call()
  ctime <- (proc.time() - ptm)[3]

  list(minimum = c(h_space, h_time), model = model, ctime = ctime)
}







