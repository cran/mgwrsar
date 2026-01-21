#' MGWR
#' to be documented
#' @usage MGWR(Y, XC, XV, ALL_X = NULL, S, H, NN, kernels, adaptive = FALSE,
#' Type = "GD", SE = FALSE, isgcv = FALSE, W = NULL, TP = NULL, Model, indexG = NULL,
#' Wd = NULL, dists = NULL,  ncore = 1, TP_estim_as_extrapol = FALSE,
#' noisland = FALSE, get_ts = FALSE, get_s = FALSE, get_Rk = get_Rk, alpha = 1)
#' @param Y  A vector
#' @param XC  A matrix with covariates with stationnary parameters
#' @param XV   A matrix with covariates with spatially varying parameters
#' @param S  A matrix with variables used in kernel
#' @param H  A vector of bandwidths
#' @param NN Number of spatial Neighbours for kernels computations.
#' @param kernels  A vector of kernel types
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel.
#' @param Type  Type of Genelarized kernel product ('GD' only
#'  spatial,'GDC' spatial + a categorical variable,
#'  'GDX' spatial + a continuous variable and other
#'   combinations like 'GDXXC','GDXCC',...)
#' @param SE  If standard error are computed, default FALSE
#' @param isgcv  leave one out cross validation, default FALSE
#' @param W  A weight matrix for spatial autocorrelation
#' @param TP  index of target points, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL.
#' @param dists  Precomputed Matrix of spatial distances, default NULL
#' @param noisland A boolean to avoid isle with no neighbours for non adaptive kernel, default FALSE
#' @param Model character containing the type of model:
#'  Possible values are "OLS", "SAR", "GWR" (default), "MGWR" ,
#'   "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv", "MGWRSAR_0_kc_kv",
#'   "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' @return a list of object for MGWRSAR wrapper
#' @noRd
MGWR <- function(Y, XC, XV, ALL_X = NULL, S, H, NN, kernels, adaptive = FALSE, Type = "GD", SE = FALSE, isgcv = FALSE, W = NULL, TP = NULL, Model, indexG = NULL, Wd = NULL, dists = NULL,ncore = 1, TP_estim_as_extrapol = FALSE, noisland = FALSE, get_ts = FALSE, get_s = FALSE, get_Rk = get_Rk, alpha = 1) {
  se = NULL
  sev = NULL
  coords = S[, 1:2]
  if(ncol(S) > 2) Z = matrix(S[, 3:ncol(S)]) else Z = NULL
  if (!is.null(XC)) XC <- as.matrix(XC)
  if (!is.null(XV)) XV <- as.matrix(XV)
  ALL_X = cbind(XC, XV)
  n <- NROW(Y)
  m <- ncol(XV)
  K <- ncol(XC)

  if (Model %in% c("MGWRSAR_0_kc_kv", "MGWRSAR_0_0_kv")) {
    PhWy = PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1, n))
    if (Model == "MGWRSAR_0_kc_kv")
      XC = cbind(XC, PhWy)
    else XC = as.matrix(PhWy, ncol = 1)
  }

  ### prep wd
  ########
  if(is.null(Wd)){
    Z = S[TP, ]
    stage1 = prep_w(H = H, kernels = kernels, Type = Type, adaptive = adaptive, dists = dists, indexG = indexG, alpha = alpha)
    indexG = stage1$indexG
    dists = stage1$dists
    Wd = stage1$Wd
    isolated_idx <- integer(0)
    if (kernels[1] != "gauss") {
      isolated_idx <- which(rowSums(Wd > 0) < ncol(XV))
      if(length(isolated_idx)>0){
        assign('isolated_idx',isolated_idx,envir=parent.frame())
      }
    }
  }


  if (Model == "MGWRSAR_1_kc_0") {

    # Variable coefficients = W %*% Y
    XV_star <- as.matrix(W %*% Y)

  } else if (Model %in% c("MGWRSAR_1_0_kv", "MGWRSAR_1_kc_kv")) {

    # Add the (W %*% Y) component to variable covariates
    XV_star <- cbind(XV, as.numeric(W %*% Y))

  } else {

    # Standard MGWR
    XV_star <- XV
  }

  # --- 1. Preparation of XV (Spatial Variables) ---
  XV_star_m <- as.matrix(XV_star)
  if (ncol(XV_star_m) > nrow(XV_star_m) && nrow(XV_star_m) < 100) {
    # Safety: If XV looks transposed (e.g., 3 rows, 1000 columns), transpose it back
    warning("Transposing XV_star_m inside R wrapper due to suspicious dimensions.")
    XV_star_m <- t(XV_star_m)
  }
  storage.mode(XV_star_m) <- "double"

  # --- 2. Preparation of XC (Constant Variables) ---
  # Ensure it's an N x K matrix (even if K=1)
  XC_m <- as.matrix(XC)
  if (nrow(XC_m) != nrow(XV_star_m)) {
    # If dimensions don't match, check if it's a vector or transposed
    if (length(XC_m) == nrow(XV_star_m)) {
      XC_m <- matrix(XC_m, ncol = 1)
    } else if (ncol(XC_m) == nrow(XV_star_m)) {
      XC_m <- t(XC_m)
    } else {
      stop("Incompatible dimensions between XC and XV.")
    }
  }
  storage.mode(XC_m) <- "double"

  # --- 3. Preparation of Y ---
  Y_v <- as.numeric(Y)
  storage.mode(Y_v) <- "double"

  # --- 4. Other matrices ---
  indexG_m <- as.matrix(indexG)
  storage.mode(indexG_m) <- "integer"

  Wd_m <- as.matrix(Wd)
  storage.mode(Wd_m) <- "double"

  TP_v <- as.integer(TP)
  storage.mode(TP_v) <- "integer"

  # --- 5. Handling Subsampling (Target Points) ---
  # If TP covers everyone, ensure alignment
  if (length(TP_v) == nrow(XV_star_m)) {
    indexG_m <- indexG_m[TP_v, , drop = FALSE]
    Wd_m     <- Wd_m[TP_v, , drop = FALSE]
  }

  # --- 6. THE CALL VIA .Call (Crucial for manual init.c) ---
  # The name here must correspond EXACTLY to the one in init.c
  res <- .Call("_mgwrsar_mgwr_beta_pivotal_qrp_mixed_cpp",
               XV_star_m, # XV
               Y_v,        # y
               XC_m,       # XC
               indexG_m,   # indexG
               Wd_m,       # Wd
               TP_v,       # TP
               get_ts,     # get_ts
               get_s,      # get_s
               get_Rk,     # get_Rk
               SE,         # get_se (Note: in your R code it was 'SE', in C++ 'get_se')
               PACKAGE = "mgwrsar" # Name of your package (very important)
  )
  res$Betac = as.numeric(res$Betac)
  names(res$Betac)=colnames(XC)

  if (SE && !isgcv) {

    output <- list(
      Betac = res$Betac,                 # Fixed coefficients
      Betav = res$Betav,                    # Local varying coef
      SEV   = res$SEV,                      # SE of varying part
      se    = as.numeric(res$se),                     # SE of fixed part
      edf   = n - res$tS - length(res$Betac),
      tS    = res$tS + length(res$Betac),
      Shat  = res$Shat                      # Hat matrix if requested
    )

  } else if (get_ts) {

    output <- list(
      Betac = res$Betac,
      Betav = res$Betav,
      SEV   = NULL,
      edf   = NULL,
      tS    = res$tS + length(res$Betac),
      Shat  = res$Shat
    )

  } else {

    output <- list(
      Betac = as.numeric(res$Betac),
      Betav = res$Betav,
      SEV   = NULL,
      edf   = NULL,
      tS    = NULL,
      Shat  = NULL
    )
  }

  return(output)
}
