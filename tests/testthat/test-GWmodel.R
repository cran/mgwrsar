# =============================================================================
# TEST GWR MGWRSAR vs GWmodels 
# =============================================================================

library(testthat)
library(mgwrsar)
source('../../tools/test_init.R')

test_that("GWR GD (fixed/adaptive) matches GWmodel within tolerance", {
  
  skip_if_not_installed("GWmodel")
  library(GWmodel)
  
  # -------------------------
  # Data (full 300 obs)
  # -------------------------
  test_data_full <- setup_test_data_full(
    n          = 300,
    lambda     = NULL,
    config_beta= 'default',
    config_snr = 0.9
  )
  
  data <- test_data_full$GD
  df   <- data$mydata
  coords <- as.matrix(data$coords)   # n x 2
  n <- nrow(coords)
  
  # GWmodel wants a SpatialPointsDataFrame
  spdf <- sp::SpatialPointsDataFrame(coords = coords, data = df)
  
  # Formula: ensure Intercept behavior aligned
  fml <- as.formula("Y ~ X1 + X2 + X3")
  
  # -------------------------
  # Bandwidths (arbitrary but stable)
  # -------------------------
  # Fixed: in distance units
  maxd <- max(as.matrix(dist(coords)))
  bw_fixed <- 0.25 * maxd  # arbitrary, but not too small
  
  # Adaptive: number of neighbors
  bw_adapt <- 80L          # arbitrary, but should be < n
  
  # -------------------------
  # Helper to align coefficient matrices
  # -------------------------
  align_coef <- function(B_mgwrsar, B_gwmodel) {
    # GWmodel returns columns like "(Intercept)", "X1", "X2", "X3"
    # mgwrsar uses "Intercept" sometimes -> harmonize names
    #colnames(B_mgwrsar)[colnames(B_mgwrsar) == "Intercept"] <- "(Intercept)"
    common <- intersect(colnames(B_mgwrsar), colnames(B_gwmodel))
    list(
      mg = B_mgwrsar[, common, drop = FALSE],
      gw = B_gwmodel[, common, drop = FALSE]
    )
  }
  
  # =============================================================================
  # Case 1: FIXED bandwidth (adaptive = FALSE)
  # =============================================================================
  
  mod_mg_fixed <- MGWRSAR(
    formula = fml,
    data    = df,
    coords  = coords,
    Model   = "GWR",
    kernels = "gauss",
    H       = bw_fixed,
    control = list(adaptive = FALSE, Type    = "GD",  NN = n)  # NN=n = full distances
  )
  
  mod_gw_fixed <- GWmodel::gwr.basic(
    formula  = fml,
    data     = spdf,
    bw       = bw_fixed,
    kernel   = "gaussian",
    adaptive = FALSE,
    longlat  = FALSE
  )
  
  # Coefs
  B_mg <- mod_mg_fixed@Betav
  B_gw <- as.matrix(mod_gw_fixed$SDF@data[, c("Intercept", "X1", "X2", "X3")])
  dimnames(B_mg)<-dimnames(B_gw)
  al <- align_coef(B_mg, B_gw)
  
  # Compare coefficients + fitted values
  expect_equal(al$mg, al$gw, tolerance = 1e-12)
  expect_equal(as.numeric(mod_mg_fixed@fit),
               as.numeric(mod_gw_fixed$SDF@data$yhat),
               tolerance = 1e-12)
  
  # =============================================================================
  # Case 2: ADAPTIVE bandwidth (adaptive = TRUE)
  # =============================================================================
  
  mod_mg_adapt <- MGWRSAR(
    formula = fml,
    data    = df,
    coords  = coords,
    Model   = "GWR",
    kernels = "gauss",
    H       = bw_adapt,
    control = list(adaptive = TRUE, Type    = "GD",NN = n) # NN >= H+2 dans tes conventions
  )
  
  mod_gw_adapt <- GWmodel::gwr.basic(
    formula  = fml,
    data     = spdf,
    bw       = bw_adapt,
    kernel   = "gaussian",
    adaptive = TRUE,
    longlat  = FALSE
  )
  
  B_mg <- mod_mg_adapt@Betav
  B_gw <- as.matrix(mod_gw_adapt$SDF@data[, c("Intercept", "X1", "X2", "X3")])
  dimnames(B_mg)<-dimnames(B_gw)
  al <- align_coef(B_mg, B_gw)
  # head(B_mg);
  # head(B_gw);
  
  expect_equal(al$mg, al$gw, tolerance = 1e-12)
  expect_equal(as.numeric(mod_mg_adapt@fit),
               as.numeric(mod_gw_adapt$SDF@data$yhat),
               tolerance = 1e-12)
})