# tools/write_initial_coef_hashes.R
# Génère les hashes "de référence" des matrices de coefficients par configuration.

# --- packages ---
stopifnot(requireNamespace("digest", quietly = TRUE))

hash_coef_matrix <- function(B, digits = 10) {
  stopifnot(is.matrix(B) || is.data.frame(B))
  B <- as.matrix(B)
  
  # Canonisation ordre
  if (!is.null(rownames(B))) B <- B[order(rownames(B)), , drop = FALSE]
  if (!is.null(colnames(B))) B <- B[, order(colnames(B)), drop = FALSE]
  
  # Canonisation type + arrondi
  storage.mode(B) <- "double"
  B <- round(B, digits = digits)
  
  # Canonisation NA / -0
  B[is.na(B)] <- NaN
  B[B == 0] <- 0
  
  digest::digest(B, algo = "sha256", serialize = TRUE)
}

# ---- fonction qui entraîne un modèle selon une config ----
fit_one <- function(cfg) {
  test_MGWR_SAR <- setup_test_data_full(n=cfg$n,lambda=cfg$lambda,config_beta = cfg$config_beta,config_snr =cfg$config_snr)
  W=NULL
  if(cfg$Type=='GD') {
    data <- test_MGWR_SAR$GD 
    Z=NULL
  } else if(cfg$Type %in% c('GDT','T')) {
    data <- test_MGWR_SAR$GDT
    Z=data$mydata$time
  }
  df   <- data$mydata
  coords<- data$coords
  fml<-data$formula
  
  if(cfg$Model %in% c('GWR','MGWR','MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')) {
    
    if(cfg$Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')) {
      W=kernel_matW(H=4,kernels='rectangle',coords=coords,NN=5,adaptive=TRUE,diagnull=TRUE)
    }
    if(cfg$Model %in% c('MGWR')) {
      fixed_vars='X1' 
      cfg$Model='GWR'
    }
    
    mymodel <- MGWRSAR(
      formula = fml,
      data    = df,
      fixed_vars=cfg$fixed_vars,
      coords  = coords,
      Model   = cfg$Model,
      kernels = cfg$kernels,
      H       = cfg$H,
      control = list(Z=Z,adaptive = cfg$adaptive, Type= cfg$Type, NN = cfg$NN,get_s=cfg$get_s,W=W, SE=cfg$SE,verbose = FALSE) 
    )
  } else { 
    ctrl <- list(adaptive = cfg$adaptive)
  
  if (!is.null(cfg$NN))
    ctrl$NN <- cfg$NN
  
  mymodel<-TDS_MGWR(
    formula      = fml,
    Model        = cfg$Model,
    data         = df,
    coords       = coords,
    kernels      = cfg$kernels,
    fixed_vars   = cfg$fixed_vars,
    control_tds  = cfg$control_tds,
    control      = list(Z=Z,adaptive = cfg$adaptive, Type= cfg$Type, NN = cfg$NN,get_s=cfg$get_s,W=W, SE=cfg$SE,verbose = FALSE) 
  )
  }
  
  mymodel
}
source('tools/configs_list_estimation.R')
source('tools/test_init.R')
if (Sys.getenv("SAVE_TESTS") == "1") res_model_coef<-list()

configs<-c(configs_gwr,configs_mgwr)
#id=names(configs)[1]
rows <- lapply(names(configs), function(id) {
  cfg <- configs[[id]]
  cat("Running:", id, "\n")
  
  model <- fit_one(cfg)
  if (Sys.getenv("SAVE_TESTS") == "1") res_model_coef[[id]]<-model
  
  B <- model@Betav
  if(cfg$SE & length(model@sev)>0) {
    sev <- model@sev
    B=B+sev
  }
  
  
  data.frame(
    config_id = id,
    n = nrow(B),
    p = ncol(B),
    digits = 10,
    hash = hash_coef_matrix(B, digits = 10),
    stringsAsFactors = FALSE
  )
})
if (Sys.getenv("SAVE_TESTS") == "1") save(res_model_coef,file='/Users/geniaux/Documents/Boulot/programme/mgwrsar2/test_version/res_model_coef.Rdata')

out <- do.call(rbind, rows)

# --- écriture du fichier de référence ---
ref_path <- file.path("tests", "testthat", "_coef_hashes.csv")
dir.create(dirname(ref_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, ref_path, row.names = FALSE)

cat("\nWrote:", ref_path, "\n")
print(out)





