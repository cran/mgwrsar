#' update_opt
#' to be documented
#' @usage update_opt(env = parent.frame())
#' @param env environment to evaluate in
#' @noRd
#' @return to be documented
update_opt <- function(env = parent.frame()){
  with(env, {
    ## choisi vks la liste de banwdith testé
    update_bandwidth_candidates()
    if (control_tds$check_pairs) vks <- vks[vks >= HKmin[k]]

    data$e0k <- data$e0 + BETA[, k] * X[, k]
    myformula_bk <- as.formula(paste0('e0k~-1+', k))


    res <- fop(foreach(v = unique(vks), .combine = "rbind", .inorder = FALSE
                       #,.export = c("myformula_bk", "data", "coords", "kernels", "Ht", "max_dist", "n", "control")
                       ), {
                          controlv <- control
                          if (v != max_dist) { #  | (v==max_dist & !is.null(Ht))
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
                           Betav <- rep(model_k@Betac, n)
                           model_k@AICc <- model_k@AIC
                         }
                         list(AICc = model_k@AICc, betav = Betav, e0 = residuals(model_k),
                              vk = v, vt = NULL, TS = as.numeric(model_k@TS), S = model_k@Shat)
                       })

    # cas général
    {

      ###
      gc(verbose = FALSE)

      # ============================================================
      # Filtrage dynamique des cas degeneres selon AICc
      # ============================================================
      AICc_vals <- as.numeric(res[, "AICc"])

      # AICc de reference = celui du modele courant, sinon mediane
      AICc_ref <- tryCatch({
        myAICc
      }, error = function(e) {
        median(AICc_vals[!is.na(AICc_vals)], na.rm = TRUE)
      })

      # Cas degeneres : AICc < 0.5 * AICc_ref (seuil arbitraire de securite)
      degenerate_idx <- which(!is.na(AICc_vals) & AICc_vals < 0.5 * AICc_ref)
      valid_idx      <- which(!is.na(AICc_vals) & !(AICc_vals < 0.5 * AICc_ref))


      # Avertissement si cas degeneres detectes
      if (length(degenerate_idx) > 0) {
        warning(sprintf(
          "[WARNING] Variable %s : %d degenerate models excluded (AICc < 0.5 * AICc_ref = %.2f) at iteration : %d",
          k, length(degenerate_idx), AICc_ref, i
        ))
      }

      # Selection du meilleur modele parmi les valides
      if (length(valid_idx) > 0) {
        mybest <- valid_idx[which.min(AICc_vals[valid_idx])]
      } else {
        warning(sprintf("[WARNING] Variable %s : all models degenerate or NA", k))
        mybest <- which.min(AICc_vals)
      }

      # ============================================================
      # Extraction du meilleur modele restant
      # ============================================================
      opt[k] <- unlist(res[mybest, "vk"])


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

    ### test en cours
    if(i==1){
      AICc_vals <- as.numeric(res[, "AICc"])
      mybest <- which.min(AICc_vals)
      AIC_ref_OLS=as.numeric(tail(res[, "AICc"],1))
      Vtemp=V
      continue=TRUE
      if(as.numeric(res[mybest, "vk"])==max_dist) while((length(Vtemp)-length(vks))>0 & continue){
        Vtemp=setdiff(Vtemp,vks)
        vks=head(Vtemp,ncore)
        res <- fop(foreach(v = unique(vks), .combine = "rbind", .inorder = FALSE
                           #,.export = c("myformula_bk", "data", "coords", "kernels", "Ht", "max_dist", "n", "control")
        ), {
          controlv <- control
          if (v != max_dist) { #  | (v==max_dist & !is.null(Ht))
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
            Betav <- rep(model_k@Betac, n)
            model_k@AICc <- model_k@AIC
          }
          list(AICc = model_k@AICc, betav = Betav, e0 = residuals(model_k),
               vk = v, vt = NULL, TS = as.numeric(model_k@TS), S = model_k@Shat)
        })
        #browser()
        if (is.null(dim(res))){
          res <- matrix(res, nrow = 1)
          colnames(res) <- c("AICc", "betav", "e0", "vk","vt","TS","S")

        }
        AICc_vals <- as.numeric(res[, "AICc"])
        mybest<-head(which(AICc_vals<AIC_ref_OLS),1)
        if(length(mybest)>0) {
          continue=FALSE
          opt[k]<-vks[mybest]

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
