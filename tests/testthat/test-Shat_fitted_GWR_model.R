# =============================================================================
# TEST gwr like
# =============================================================================

library(testthat)
library(mgwrsar)

source('../../tools/configs_list_estimation.R')
source('../../tools/test_init.R')
source('../../tools/check_hash_against_registry.R')

## vérifier pour tout modèle la validité de S, de fit, de résiduals sur un DGP MGTWR_SAR n=600

test_that("GWR like model : test S fit residuals coherence", {
  configs<-configs_gwr
  # -------------------------
  # Build data ONCE (all your configs here use the same DGP)
  # -------------------------
  # Sanity check: all configs share these DGP params
  cfg0 <- configs[[1]]
  stopifnot(
    all(vapply(configs, function(cfg) cfg$n, numeric(1)) == cfg0$n),
    all(vapply(configs, function(cfg) cfg$lambda, numeric(1)) == cfg0$lambda),
    all(vapply(configs, function(cfg) cfg$config_beta, character(1)) == cfg0$config_beta),
    all(vapply(configs, function(cfg) cfg$config_snr, numeric(1)) == cfg0$config_snr),
    all(vapply(configs, function(cfg) cfg$Type, character(1)) == "GD")
  )
  
  test_MGWR_SAR <- setup_test_data_full(
    n          = cfg0$n,
    lambda     = cfg0$lambda,
    config_beta= cfg0$config_beta,
    config_snr = cfg0$config_snr
  )
  
  data   <- test_MGWR_SAR$GD
  df     <- data$mydata
  coords <- as.matrix(data$coords)
  fml    <- data$formula
  
  # Same W as your original test
  W <- kernel_matW(
    H = 4, kernels = "rectangle",
    coords = coords, NN = 5, adaptive = TRUE,
    diagnull = TRUE
  )
  
  # -------------------------
  # Helper to fit one config
  # -------------------------
  fit_one <- function(cfg, get_s) {
    MGWRSAR(
      formula    = fml,
      data       = df,
      fixed_vars = cfg$fixed_vars,
      coords     = coords,
      Model      = cfg$Model,
      kernels    = cfg$kernels,
      H          = cfg$H,
      control    = list(
        adaptive = cfg$adaptive,
        Type     = cfg$Type,
        NN       = cfg$NN,
        get_s    = get_s,
        W        = W
      )
    )
  }
  
  # -------------------------
  # Run all configs
  # -------------------------
  for (id in names(configs)) {
    cfg <- configs[[id]]
    #cat(id,'\n')
    
    # two fits: without S and with S
    mymodel   <- fit_one(cfg, get_s = FALSE)
    mymodel_S <- fit_one(cfg, get_s = TRUE)
  
    
    # When Shat is returned, trace(Shat) must match tS from the non-S fit
    expect_equal(sum(diag(mymodel_S@Shat)), mymodel@tS, tolerance = 1e-12,
                 info = paste("trace(Shat) mismatch for", id))
    
    # Fitted/residuals must match between get_s=F and get_s=T
    expect_equal(fitted(mymodel_S), fitted(mymodel), tolerance = 1e-12,
                 info = paste("fitted mismatch for", id))
    expect_equal(residuals(mymodel_S), residuals(mymodel), tolerance = 1e-12,
                 info = paste("residuals mismatch for", id))
    expect_equal(df$Y, fitted(mymodel_S) + residuals(mymodel_S), tolerance = 1e-12,
                 info = paste("Y != fitted + residuals for", id))
    
    # Extra check for pure GWR (i.e., no SAR operator complications)
    # Note: your original code only did this for GWR/MGWR-like (Model == "GWR")
    if (identical(cfg$Model, "GWR")) {
      expect_equal(as.numeric(mymodel_S@Shat %*% df$Y), fitted(mymodel), tolerance = 1e-12,
                   info = paste("Shat %*% Y mismatch for", id))
    }
    
  
    B <- mymodel@Betav
    if(cfg$SE & length(mymodel@sev)>0) {
      sev <- mymodel@sev
      B=B+sev
    }
    
    row=data.frame(
      config_id = id,
      n = nrow(B),
      p = ncol(B),
      digits = 10,
      hash = hash_coef_matrix(B, digits = 10),
      stringsAsFactors = FALSE
    )
    check_hash_against_registry_coef(row,file='_coef_hashes.csv')
  }
  
})
