gwr_beta_pivotal_qrp_full <- function(
    Y, XV, ALL_X, TP, indexG, Wd, NN,
    W = NULL,
    isgcv = FALSE, SE = FALSE,
    kernels = NULL, H = NULL, adaptive = NULL,
     ncore = 1,
    TP_estim_as_extrapol = FALSE,
    get_ts = FALSE, get_s = FALSE, get_Rk = FALSE,
    isolated_idx = NULL
) {
  if (get_s) get_ts <- TRUE
  if (isgcv) Wd[,1] <- 0

  n <- length(Y)
  ntp <- length(TP)
  m <- if (!is.null(XV)) ncol(XV) else 0
  namesXV <- colnames(XV)

  if (!is.null(isolated_idx)) TP <- isolated_idx

  if (TP_estim_as_extrapol) {
    SE <- FALSE
    isgcv <- FALSE
  }

  # --- Initialization of containers
  if (get_s) Shat <- matrix(0, n, n) else Shat <- NULL
  if (get_s | SE) SEV <- matrix(0, n, ncol(XV)) else SEV <- NULL
  TS <- rep(0, n)
  tS <- 0
  edf <- n
  Rk <- NULL   # <-- important to avoid "object 'Rk' not found" error

  # --- Case W (SAR)
  if (!is.null(W)) {
    PhWy <- PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1, n))
    XV <- cbind(XV, PhWy)
  }

  # --- Automatic choice: univariate or multivariate
  if (ncol(XV) == 1L) {
    res_cpp <- gwr_beta_univar_cpp(
      y = Y,
      x = XV[, 1],
      XV = XV,
      indexG = indexG,
      Wd = Wd,
      TP = TP,
      get_ts = get_ts,
      get_s = get_s
    )
  } else {
    res_cpp <- gwr_beta_pivotal_qrp_cpp(
      X = XV,
      y = Y,
      XV = XV,
      indexG = indexG,
      Wd = Wd,
      TP = TP,
      get_ts = get_ts,
      get_s = get_s,
      get_Rk = get_Rk,
      get_se = SE
    )
  }

  # --- Result extraction
  Betav <- res_cpp$Betav
  if (anyNA(Betav) || any(rowSums(abs(Betav)) == 0)) {

    # Global OLS estimation (robust fallback)
    beta_ols <- tryCatch(
      coef(lm.fit(x = ALL_X, y = Y)),
      error = function(e) rep(NA_real_, ncol(ALL_X))
    )

    if (!anyNA(beta_ols)) {
      # Detect rows to replace
      replace_idx <- which(
        apply(Betav, 1, function(x) all(is.na(x)) || sum(abs(x)) < .Machine$double.eps)
      )
      if (length(replace_idx) > 0) {
        Betav[replace_idx, ] <- matrix(
          beta_ols,
          nrow = length(replace_idx),
          ncol = length(beta_ols),
          byrow = TRUE
        )
      }
    }
  }

  # Update main object
  res_cpp$Betav <- Betav
  if (get_ts) TS[TP] <- res_cpp$TS
  if (get_s)  Shat   <- res_cpp$Shat

  # --- Secure conversion of Rk
  if (get_Rk && "Rk" %in% names(res_cpp)) {
    Rk <- res_cpp$Rk
    if (!is.list(Rk)) {
      # convert array [nTP, n, p] -> list of matrices
      p <- dim(Rk)[3]
      Rk_list <- vector("list", p)
      for (k in seq_len(p)) Rk_list[[k]] <- Rk[, , k]
      Rk <- Rk_list
      names(Rk) = colnames(XV)
    }
  }

  if (get_ts) tS <- sum(TS, na.rm = TRUE)
  edf <- n - tS

  # --- Local SE (optional)
  if (SE && !isgcv) {
    SEV = res_cpp$SEV
  }

  if (!is.null(W))
    colnames(Betav) <- c(namesXV, "lambda")
  else
    colnames(Betav) <- namesXV

  # --- Final result construction
  out <- list(Betav = Betav)

  if (get_s | get_ts | SE) {
    out$SEV  <- SEV
    out$edf  <- edf
    out$tS   <- tS
    out$Shat <- Shat
    out$TS   <- TS
  }

  if (get_Rk) out$Rk <- Rk

  return(out)
}
