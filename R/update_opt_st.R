#' Update Optimal Bandwidth for a Single Variable (Spatio-Temporal)
#'
#' Evaluates candidate spatial and temporal bandwidth pairs for a single
#' variable in the backfitting procedure and selects the pair minimizing AICc.
#'
#' @usage update_opt_st(env = parent.frame())
#' @param env Environment to evaluate in (typically the parent backfitting frame).
#' @noRd
#' @return Modifies variables in the parent environment (opt, opt_t, e0, betav, etc.).
update_opt_st <- function(env = parent.frame()) {
  with(env, {
    
    ncore_use <- control_tds$ncore
   # if (Sys.info()[["sysname"]] == "Darwin") {
  #    ncore_use <- 1L
   # }

    update_bandwidth_candidates()

    data$e0k <- data$e0 + BETA[, k] * X[, k]
    myformula_bk <- as.formula(paste0("e0k ~ -1 + ", k))

    vkst <- unique(expand.grid(unique(vks), unique(vkt)))

    if(control_tds$check_pairs){
      hmin = HKMIN[[k]]
      # Protection: match can return NA if no correspondence, so we filter
      idx <- match(vkst[, 1], hmin[, 1])
      vkst <- vkst[which(vkst[, 2] >= hmin[idx, 2]), , drop = FALSE]
    }

    # controlvd <- modifyList(control, list(
    #   dists = NULL, indexG = NULL, Type = "GD",
    #   adaptive = control$adaptive[1]
    # ))
    #
    # controlvt <- modifyList(control, list(
    #   dists = NULL, indexG = NULL, Type = "T",
    #   adaptive = FALSE
    # ))

    .eval_bw_st <- function(j) {
      controlv <- control
      v  <- vkst[j, 1]
      vt <- vkst[j, 2]
      model_k <- NULL
      if (v < max_dist & vt < max_dist_t) {
        model_k <- MGWRSAR(
          formula = myformula_bk, data = data, coords = coords,
          fixed_vars = NULL, kernels = kernels, H = c(v, vt),
          Model = "GWR", control = controlv
        )
      } else if (v == max_dist & vt == max_dist_t) {
        model_k <- MGWRSAR(
          formula = myformula_bk, data = data, coords = coords,
          fixed_vars = NULL, kernels = kernels, H = NULL,
          Model = "OLS", control = controlv
        )
        model_k@sev=as.matrix(rep(model_k@se,n),ncol=1)
        model_k@Betav <- as.matrix(rep(model_k@Betac,n),ncol=1)
        model_k@AICc <- model_k@AIC
      } else if (v < max_dist & vt == max_dist_t) {
        model_k <- MGWRSAR(
          formula = myformula_bk, data = data, coords = coords,
          fixed_vars = NULL, kernels = kernels[1], H = c(v, NULL),
          Model = "GWR", control = controlvd
        )
      } else if (v == max_dist & vt < max_dist_t) {
        model_k <- MGWRSAR(
          formula = myformula_bk, data = data, coords = as.matrix(control$Z, ncol = 1),
          fixed_vars = NULL, kernels = kernels[2], H = vt,
          Model = "GWR", control = controlvt
        )
      }
      list(
        AICc  = model_k@AICc,
        betav = model_k@Betav,
        e0    = residuals(model_k),
        vk    = v,
        vt    = vt,
        TS    = as.numeric(model_k@TS),
        S     = model_k@Shat,
        sev   = model_k@sev,
        edf   = model_k@edf
      )
    }

    if (ncore_use > 1) {
      res <- do.call(rbind, parallel::mclapply(1:nrow(vkst), .eval_bw_st, mc.cores = ncore_use))
    } else {
      res <- do.call(rbind, lapply(1:nrow(vkst), .eval_bw_st))
    }

    AICc_vals <- as.numeric(res[, "AICc"])

    AICc_ref <- tryCatch({
      myAICc
    }, error = function(e) {
      median(AICc_vals[!is.na(AICc_vals)], na.rm = TRUE)
    })

    degenerate_idx <- which(!is.na(AICc_vals) & AICc_vals < 0.5 * AICc_ref)
    valid_idx <- which(!is.na(AICc_vals) & !(AICc_vals < 0.5 * AICc_ref))

    if (length(degenerate_idx) > 0) {
      warning(sprintf(
        "Variable %s : %d degenerate models excluded (AICc < 0.5 * AICc_ref = %.2f) at iteration : %d",
        k, length(degenerate_idx), AICc_ref, i
      ))
    }

    if (length(valid_idx) > 0) {
      mybest <- valid_idx[which.min(AICc_vals[valid_idx])]
    } else {
      warning(sprintf("Variable %s : all models degenerate or NA", k))
      mybest <- which.min(AICc_vals)
    }

    opt[k]   <- unlist(res[mybest, "vk"])
    opt_t[k] <- unlist(res[mybest, "vt"])

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
  })
}

