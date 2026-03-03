#' Update Optimal Bandwidth for a Single Variable (Spatial)
#'
#' Evaluates candidate bandwidths for a single variable in the backfitting
#' procedure and selects the one minimizing AICc with degenerate case filtering.
#'
#' @usage update_opt(env = parent.frame())
#' @param env Environment to evaluate in (typically the parent backfitting frame).
#' @noRd
#' @return Modifies variables in the parent environment (opt, e0, betav, etc.).
update_opt <- function(env = parent.frame()){
  with(env, {
    
    # ------------------------------------------------------------
    # FORCE SEQUENTIAL ON macOS (avoid fork deadlock)
    # ------------------------------------------------------------
    ncore_use <- control_tds$ncore
  #  if (Sys.info()[["sysname"]] == "Darwin") {
   #   ncore_use <- 1L
   # }
    
    ## select vks the list of bandwidth candidates to test
    update_bandwidth_candidates()
    if (control_tds$check_pairs) vks <- vks[vks >= HKmin[k]]
    
    data$e0k <- data$e0 + BETA[, k] * X[, k]
    myformula_bk <- as.formula(paste0('e0k~-1+', k))
    
    
    .eval_bw <- function(v) {
      controlv <- control
      if (v != max_dist) {
        if (control$adaptive[1]) {
          if (kernels[1] != 'gauss') {
            NNN <- v + 2
            controlv$NN <- min(NNN, control$NN)
            controlv$indexG <- control$indexG[, 1:controlv$NN]
            controlv$dists[['dist_s']] <- control$dists[['dist_s']][, 1:controlv$NN]
          }
        }
        model_k <- MGWRSAR(formula = myformula_bk, data = data, coords = coords,
                           fixed_vars = NULL, kernels = kernels, H = c(v, NULL),
                           Model = 'GWR', control = controlv)
        Betav <- model_k@Betav
      } else {
        model_k <- MGWRSAR(formula = myformula_bk, data = data, coords = coords,
                           fixed_vars = NULL, kernels = kernels, H = c(v, NULL),
                           Model = 'OLS', control = controlv)
        model_k@sev = as.matrix(rep(model_k@se,n),ncol=1)
        Betav <- rep(model_k@Betac,n)
        model_k@AICc <- model_k@AIC
      }
      list(AICc = model_k@AICc, betav = Betav, e0 = residuals(model_k),
           vk = v, vt = NULL, TS = as.numeric(model_k@TS),
           S = model_k@Shat, sev=model_k@sev, edf=model_k@edf)
    }
    
    # ------------------------------------------------------------
    # Parallel only if not macOS
    # ------------------------------------------------------------
    if (ncore_use > 1) {
      res <- do.call(rbind,
                     parallel::mclapply(unique(vks), .eval_bw,
                                        mc.cores = ncore_use,
                                        mc.preschedule = FALSE))
    } else {
      res <- do.call(rbind, lapply(unique(vks), .eval_bw))
    }
    
    # ============================================================
    # general case
    # ============================================================
    {
      AICc_vals <- as.numeric(res[, "AICc"])
      
      AICc_ref <- tryCatch({
        myAICc
      }, error = function(e) {
        median(AICc_vals[!is.na(AICc_vals)], na.rm = TRUE)
      })
      
      degenerate_idx <- which(!is.na(AICc_vals) &
                                AICc_vals < 0.5 * AICc_ref)
      valid_idx <- which(!is.na(AICc_vals) &
                           !(AICc_vals < 0.5 * AICc_ref))
      
      if (length(degenerate_idx) > 0) {
        warning(sprintf(
          "[WARNING] Variable %s : %d degenerate models excluded (AICc < 0.5 * AICc_ref = %.2f) at iteration : %d",
          k, length(degenerate_idx), AICc_ref, i
        ))
      }
      
      if (length(valid_idx) > 0) {
        mybest <- valid_idx[which.min(AICc_vals[valid_idx])]
      } else {
        warning(sprintf("[WARNING] Variable %s : all models degenerate or NA", k))
        mybest <- which.min(AICc_vals)
      }
      
      opt[k] <- unlist(res[mybest, "vk"])
      
      if (any(stable > 0) & !spacestable)
        spacestable <- TRUE
      
      e0 <- unlist(res[mybest, "e0"])
      betav <- unlist(res[mybest, "betav"])
      isol <- is.na(betav)
      e0[isol] <- data$e0[isol]
      
      if(control$SE) {
        SEV[,k]= unlist(res[mybest, "sev"])
        EDF[k]= unlist(res[mybest, "edf"])
      }
      
      if (get_AIC) {
        TSik[!isol, k] <- unlist(res[mybest, "TS"])[!isol]
        Sk <- unlist(res[mybest, "S"][[1]])
      }
    }
    
    # ============================================================
    # work in progress
    # ============================================================
    if(i==1){
      AICc_vals <- as.numeric(res[, "AICc"])
      mybest <- which.min(AICc_vals)
      AIC_ref_OLS = as.numeric(tail(res[, "AICc"],1))
      Vtemp = V
      continue = TRUE
      
      if(as.numeric(res[mybest, "vk"]) == max_dist)
        while((length(Vtemp)-length(vks))>0 & continue){
          
          Vtemp=setdiff(Vtemp,vks)
          
          # FIX: use ncore_use
          vks=head(Vtemp,ncore_use)
          
          .eval_bw2 <- function(v) {
            controlv <- control
            if (v != max_dist) {
              if (control$adaptive[1]) {
                if (kernels[1] != 'gauss') {
                  NNN <- v + 2
                  controlv$NN <- min(NNN, control$NN)
                  controlv$indexG <- control$indexG[, 1:controlv$NN]
                  controlv$dists[['dist_s']] <- control$dists[['dist_s']][, 1:controlv$NN]
                }
              }
              model_k <- MGWRSAR(formula = myformula_bk, data = data,
                                 coords = coords,
                                 fixed_vars = NULL,
                                 kernels = kernels,
                                 H = c(v, NULL),
                                 Model = 'GWR',
                                 control = controlv)
              Betav <- model_k@Betav
            } else {
              model_k <- MGWRSAR(formula = myformula_bk, data = data,
                                 coords = coords,
                                 fixed_vars = NULL,
                                 kernels = kernels,
                                 H = c(v, NULL),
                                 Model = 'OLS',
                                 control = controlv)
              Betav <- rep(model_k@Betac, n)
              model_k@AICc <- model_k@AIC
            }
            list(AICc = model_k@AICc, betav = Betav,
                 e0 = residuals(model_k),
                 vk = v, vt = NULL,
                 TS = as.numeric(model_k@TS),
                 S = model_k@Shat)
          }
          
          if (ncore_use > 1) {
            res <- do.call(rbind,
                           parallel::mclapply(unique(vks),
                                              .eval_bw2,
                                              mc.cores = ncore_use,
                                              mc.preschedule = FALSE))
          } else {
            res <- do.call(rbind,
                           lapply(unique(vks), .eval_bw2))
          }
          
          if (is.null(dim(res))){
            res <- matrix(res, nrow = 1)
            colnames(res) <- c("AICc","betav","e0","vk","vt","TS","S")
          }
          
          AICc_vals <- as.numeric(res[, "AICc"])
          mybest <- head(which(AICc_vals < AIC_ref_OLS),1)
          
          if(length(mybest)>0) {
            continue=FALSE
            opt[k]=vks[mybest]
            
            if (any(stable > 0) & !spacestable)
              spacestable <- TRUE
            
            e0 <- unlist(res[mybest, "e0"])
            betav <- unlist(res[mybest, "betav"])
            isol <- is.na(betav)
            e0[isol] <- data$e0[isol]
            
            if (get_AIC) {
              TSik[!isol, k] <- unlist(res[mybest, "TS"])[!isol]
              Sk <- unlist(res[mybest, "S"][[1]])
            }
          }
        }
    }
    
  })
}
