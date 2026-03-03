check_hash_against_registry <- function(row, file = "results_registry.csv") {
  
  if (!file.exists(file)) {
    stop("Registry file does not exist.")
  }
  
  registry <- read.csv(file, stringsAsFactors = FALSE)
  
  # Identifier la ligne correspondant à la configuration
  match_idx <- with(registry,
                    config_id   == row$config_id   &
                      cfg_name    == row$cfg_name    &
                      Model       == row$Model       &
                      kernels     == row$kernels     &
                      adaptive    == row$adaptive    &
                      criterion   == row$criterion   &
                      NNtrain     == row$NNtrain     &
                      n           == row$n           &
                      hs_range    == row$hs_range    &
                      method_pred == row$method_pred &
                      digits      == row$digits
  )
  
  if (!any(match_idx)) {
    stop("No matching configuration found in registry.")
  }
  
  if (sum(match_idx) > 1) {
    stop("Duplicate configuration rows found in registry.")
  }
  
  stored_hash <- registry$hash_pred[match_idx]
  current_hash <- row$hash_pred
  
  norm_hash <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\r", "", x, fixed = TRUE)
    x
  }
  
  stored_hash  <- norm_hash(registry$hash_pred[match_idx])[1]
  current_hash <- norm_hash(row$hash_pred)[1]
  
  
  if (!identical(stored_hash, current_hash)) {
    cat("Stored: ", stored_hash, "\n")
    cat("Current:", current_hash, "\n")
    cat("Stored raw: "); print(charToRaw(stored_hash))
    cat("Curr   raw: "); print(charToRaw(current_hash))
    stop("Hash mismatch detected!")
  }
  
  message("Hash matches registry ✓")
  return(TRUE)
}

check_hash_against_registry_coef <- function(row, file = "results_registry.csv") {
  
  if (!file.exists(file)) {
    stop("Registry file does not exist.")
  }
  
  registry <- read.csv(file, stringsAsFactors = FALSE)
  
  # Identifier la ligne correspondant à la configuration
  match_idx <- with(registry,
                    config_id == row$config_id &
                    n == row$n &
                    p == row$p &
                    digits == row$digits
                    
  )
  
  if (!any(match_idx)) {
    stop("No matching configuration found in registry.")
  }
  
  if (sum(match_idx) > 1) {
    stop("Duplicate configuration rows found in registry.")
  }
  
  stored_hash <- registry$hash[match_idx]
  current_hash <- row$hash
  
  norm_hash <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\r", "", x, fixed = TRUE)
    x
  }
  
  stored_hash  <- norm_hash(registry$hash_pred[match_idx])[1]
  current_hash <- norm_hash(row$hash_pred)[1]
  
  
  if (!identical(stored_hash, current_hash)) {
    cat("Stored: ", stored_hash, "\n")
    cat("Current:", current_hash, "\n")
    cat("Stored raw: "); print(charToRaw(stored_hash))
    cat("Curr   raw: "); print(charToRaw(current_hash))
    stop("Hash mismatch detected!")
  }
  
  message("Hash matches registry ✓")
  return(TRUE)
}