# =============================================================================
# TEST mgwr like
# =============================================================================

library(testthat)
library(mgwrsar)
source('../../tools/configs_list_estimation.R')
source('../../tools/test_init.R')
source('../../tools/check_hash_against_registry.R')
#source('tools/configs_list_estimation.R')
#source('tools/test_init.R')


## vérifier pour tout modèle la validité de S, de fit, de résiduals sur un DGP MGTWR_SAR n=600
test_that("MGWR like model : test S fit residuals coherence", {
  
  if (Sys.getenv("RUN_LONG_TESTS") == "1") configs_mgwr<-configs_mgwr[c(7,10)] else  if (Sys.getenv("RUN_LONG_TESTS") == "2") configs_mgwr<-configs_mgwr[c(7,10,13)] else  if (Sys.getenv("RUN_LONG_TESTS") == "3") configs_mgwr<-configs_mgwr[1:length(configs_mgwr)] else skip("MGWR and MGWR tests disabled (set RUN_LONG_TESTS='1' for single MGWR,'2' for single MGWR + MGTWR), '3' for FULL test ")
  
  # -------------------------
  # Build data once
  # -------------------------
  
  
  # -------------------------
  # Helper fit function
  # -------------------------
  fit_one <- function(cfg) {
    
    if(cfg$Type=='GDT') Z=df$time else Z=NULL
    ctrl <- list(Z=Z,SE=cfg$SE,Type=cfg$Type,NN=cfg$NN,adaptive=cfg$adaptive)
    
    if (!is.null(cfg$NN))
      ctrl$NN <- cfg$NN
    
    TDS_MGWR(
      formula      = cfg$fml,
      Model        = cfg$Model,
      data         = df,
      coords       = coords,
      kernels      = cfg$kernels,
      fixed_vars   = cfg$fixed_vars,
      control_tds  = cfg$control_tds,
      control      = ctrl
    )
  }
  
  
  # -------------------------
  # Loop over configs
  # -------------------------
  for (id in names(configs_mgwr)) {
    
    cfg <- configs_mgwr[[id]]
    #cat(id,'\n')
    

    test_MGWR_SAR <- setup_test_data_full(
      n          = cfg$n,
      lambda     = cfg$lambda,
      config_beta= cfg$config_beta,
      config_snr = cfg$config_snr
    )
    
    data   <- test_MGWR_SAR$GD
    df     <- data$mydata
    coords <- as.matrix(data$coords)
    fml<-cfg$fml    <- data$formula
    
    model <- fit_one(cfg)
 
    # z=apply(model@Betav-df[,5:8] %>% as.data.frame(),2,rmse);z;mean(z);
    
    # ---- tds_mgwr coherence ----
    if (cfg$Model == "tds_mgwr" && isTRUE(cfg$control_tds$get_AIC)) {
      
      cfg_noS <- cfg
      cfg_noS$control_tds$get_AIC <- FALSE
      
      model_noS <- fit_one(cfg_noS)
      
      expect_equal(as.numeric(model@Shat %*% df$Y),
                   as.numeric(fitted(model_noS)),
                   tolerance = 1e-12,
                   info = id)
      
      expect_equal(sum(diag(model@Shat)),
                   model@tS,
                   tolerance = 1e-12,
                   info = id)
      
      expect_equal(fitted(model),
                   fitted(model_noS),
                   tolerance = 1e-12,
                   info = id)
      
      expect_equal(residuals(model),
                   residuals(model_noS),
                   tolerance = 1e-12,
                   info = id)
    }
    
    # ---- atds_mgwr coherence ----
    if (cfg$Model == "atds_mgwr") {
      
      expect_equal(as.numeric(model@Shat %*% df$Y),
                   fitted(model),
                   tolerance = 1e-4,
                   info = id)
      
      expect_equal(sum(diag(model@Shat)),
                   model@tS,
                   tolerance = 1e-12,
                   info = id)
      
      expect_equal(as.numeric(df$Y),
                   as.numeric(fitted(model) + residuals(model)),
                   tolerance = 1e-12,
                   info = id)
    }
    
    B <- model@Betav
    if(cfg$SE & length(model@sev)>0) {
      sev <- model@sev
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

