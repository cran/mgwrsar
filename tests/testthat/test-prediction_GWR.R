library(testthat)
library(mgwrsar)

source('../../tools/configs_list_prediction_GWR.R')
source('../../tools/check_hash_against_registry.R')
source('../../tools/test_init.R')

test_that("prediction test GWR ", {
  
  if (Sys.getenv("SAVE_TESTS") == "1") res_model_GWR<-res_pred_GWR<-list()
  
  stopifnot(requireNamespace("digest", quietly = TRUE))
  
  stopifnot(exists("configs_pred_GWR"))
  configs <- configs_pred_GWR
  
  # -------------------------
  # Hash helper
  # -------------------------
  hash_pred <- function(y, digits = 12) {
    y <- as.numeric(y)
    storage.mode(y) <- "double"
    y <- round(y, digits = digits)
    y[is.na(y)] <- NaN
    y[y == 0] <- 0
    digest::digest(y, algo = "sha256", serialize = TRUE)
  }
  
  # -------------------------
  # CSV management
  # -------------------------
  ref_path <- file.path("tests", "testthat", "_pred_hashes.csv")
  dir.create(dirname(ref_path), recursive = TRUE, showWarnings = FALSE)
  
  write_or_merge_row <- function(row_df) {
    if (file.exists(ref_path)) {
      old <- read.csv(ref_path, stringsAsFactors = FALSE)
      out <- rbind(old, row_df)
      out <- out[!duplicated(out$config_id, fromLast = TRUE), , drop = FALSE]
      out <- out[order(out$config_id), , drop = FALSE]
    } else {
      out <- row_df[order(row_df$config_id), , drop = FALSE]
    }
    write.csv(out, ref_path, row.names = FALSE)
  }
  
  rmse <- function(observed, predicted) sqrt(mean((observed - predicted)^2))
  
  # ============================================================
  # Loop over configs
  # ============================================================
  for (cfg_name in names(configs)) {
    
    cfg <- configs[[cfg_name]]
    #cat('cfg_name',cfg_name,'\n')
    # ---- data generation ----
    test_MGWR_SAR <- setup_test_data_full(
      n = cfg$n,
      lambda = cfg$lambda,
      config_beta = cfg$config_beta,
      config_snr = cfg$config_snr
    )
    
    if (is.null(cfg$Type) || cfg$Type == "GD") {
      data <- test_MGWR_SAR$GD
      Z <- NULL
    } else if (cfg$Type %in% c("GDT", "T")) {
      data <- test_MGWR_SAR$GDT
      Z <- data$mydata$time
    } else {
      stop("Unknown cfg$Type: ", cfg$Type)
    }
    
    df     <- data$mydata
    coords <- as.matrix(data$coords)
    n      <- nrow(coords)
    fml    <- data$formula
    
    W_in<-NULL
    # ---- fixed test split (kt removed from config) ----
    kt <- 20
    id_test  <- find_TP(fml, df, coords, kt = kt)
    id_train <- setdiff(seq_len(n), id_test)
    
    NN_train <- if (!is.null(cfg$NN_train)) cfg$NN_train else length(id_train)
    
    newdata        <- df[id_test, , drop = FALSE]
    newdata_coords <- coords[id_test, , drop = FALSE]
    
    # ---- estimation ----
    kernels_est   <- cfg$kernels_est
    Model_est     <- cfg$Model_est
    adaptive_est  <- cfg$adaptive_est
    
    if(Model_est %in% c('GWR','MGWR','MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')){
      criterion_est <- if (!is.null(cfg$criterion_est)) cfg$criterion_est else "AICc"
      
      if(Model_est %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')){
        W_in_out=kernel_matW(H=4,kernels='rectangle',coords=rbind(coords[id_train,],coords[id_test,]),NN=4,adaptive=TRUE,diagnull=TRUE)
        W_in=W_in_out[id_train,id_train]
        W_in=mgwrsar::normW(W_in)
        Type=cfg$Type
      }
      
      hs_range_est <- cfg$hs_range
      if (is.na(hs_range_est[2]) | hs_range_est[2]>NN_train) hs_range_est[2] <- NN_train
      
      sb <- search_bandwidths(
        formula = fml,
        data    = df[id_train,],
        coords  = coords[id_train,],
        kernels = kernels_est,
        Model   = Model_est,
        control = list(
          Z = Z,
          adaptive  = adaptive_est,
          criterion = criterion_est,
          verbose   = FALSE,
          NN        = NN_train,
          W         = W_in
        ),
        hs_range = hs_range_est
      )
      
      mymodel <- sb$best_model
    } else if(Model_est %in% c("tds_mgwr","atds_mgwr")) {
      mymodel <-TDS_MGWR(
        formula      = fml,
        Model        = Model_est,
        data         = df[id_train,],
        coords       = coords[id_train,],
        kernels      = kernels_est,
        control_tds  = list(nns=20, get_AIC=FALSE, verbose=FALSE, ncore=8),
        control      = list(adaptive  = adaptive_est,verbose= FALSE,NN= NN_train)
      )
    }
    methods<- cfg$methods
    #methods <- if (!is.null(cfg$methods)) cfg$methods else c("shepard", "tWtp_model", "model", "TP")
    if (Sys.getenv("SAVE_TESTS") == "1") res_model_GWR[[cfg_name]]<-mymodel
    
    # ---- ID builder (kt + Hbest removed) ----
    build_id <- function(method_pred) {
      paste0(
        "pred",
        "::cfg=", cfg_name,
        "::Model=", Model_est,
        "::kern=", paste(kernels_est, collapse = "+"),
        "::adaptive=", adaptive_est,
        "::criterion=", criterion_est,
        "::NNtrain=", NN_train,
        "::n=", n,
        "::hsrange=", paste(hs_range_est, collapse = "-"),
        "::method=", method_pred
      )
    }
    
    
    # ---- prediction loop ----
    for (method_pred in methods) {
      if(Model_est %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')){
        Y_pred <- predict(
          mymodel,
          newdata = newdata,
          newdata_coords = newdata_coords,
          W=W_in_out,
          method_pred=method_pred,
          type= Type
        )
      } else {
        Y_pred <- predict(
          mymodel,
          newdata = newdata,
          newdata_coords = newdata_coords,
          method_pred = method_pred
        )
      }
      
      rmse_pred <- rmse(df$Y[id_test], Y_pred)
      expect_true(rmse_pred / mymodel@RMSE < 4)
      if (Sys.getenv("SAVE_TESTS") == "1")  res_model_GWR[[paste0(cfg_name,'_',method_pred)]]<-Y_pred
      
      row <- data.frame(
        config_id   = build_id(method_pred),
        cfg_name    = cfg_name,
        Model       = Model_est,
        kernels     = paste(kernels_est, collapse = "+"),
        adaptive    = adaptive_est,
        criterion   = criterion_est,
        NNtrain     = NN_train,
        n           = n,
        hs_range    = paste(hs_range_est, collapse = "-"),
        method_pred = method_pred,
        digits      = 12,
        hash_pred   = hash_pred(Y_pred, digits = 12),
        stringsAsFactors = FALSE
      )
      
      check_hash_against_registry(row,file='_pred_hashes.csv')
      
    }
  }
  if (Sys.getenv("SAVE_TESTS") == "1")  save(res_model_GWR,res_pred_GWR,file='/Users/geniaux/Documents/Boulot/programme/mgwrsar2/test_version/res_GWR.Rdata')
  
})
