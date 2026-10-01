update_bandwidth_candidates <- function(env = parent.frame()) {
  with(env, {
    # Safety: default values
    if (!exists("V5")) V5 <- NULL
    if (!exists("V5t")) V5t <- NULL

    # -----------------------------------------------------------------
    #  Spatial part: compute vks
    # -----------------------------------------------------------------
    if (pinned_s[k]) {
      # Bandwidth pinned by control_tds$H: a single candidate. up/down/ddown are
      # left untouched so that a pinned coefficient does not feed the ddown of
      # the free ones.
      vks <- opt[k]
    } else {
      # V may stop below max_dist (e.g. control_tds$first_nn < n in adaptive mode),
      # in which case there is no candidate above opt[k] and tail() returns a
      # zero-length value: keep the current position rather than assigning it.
      if (opt[k] < max_dist) {
        up_cand <- V[V > opt[k]]
        up[k] <- if (length(up_cand)) tail(up_cand, 1) else opt[k]
      } else {
        up[k] <- max_dist
      }

      if (opt[k] > min_dist) {
        if (sum(V < opt[k]) > 0)
          down[k] <- head(V[V < opt[k]], 1)
        else
          down[k] <- min_dist
        ddown[k] <- max(min(down[down > min_dist], ddown[k], na.rm = TRUE),
                        min_dist, na.rm = TRUE)
      } else {
        ddown[k] <- down[k] <- min_dist
      }

      # --- Add extended candidates V5 ---
      if (!is.null(V5)) {
        v5_up <- V5[V5 > up[k]]
        if (length(v5_up) > 0) {
          up5 <- min(v5_up)
          extra_v <- V5[V5 >= up5]
          vks <- sort(unique(c(ddown[k], down[k], opt[k], up[k], extra_v)))
        } else {
          vks <- sort(unique(c(ddown[k], down[k], opt[k], up[k])))
        }
      } else {
        vks <- sort(unique(c(ddown[k], down[k], opt[k], up[k])), decreasing = TRUE)
      }
    }

    # -----------------------------------------------------------------
    #  Temporal part: compute vkt (if opt_t and Vt exist)
    # -----------------------------------------------------------------
    if (exists("opt_t") && !is.null(opt_t) && exists("Vt") && !is.null(Vt)) {
      if (pinned_t[k]) {
        # Bandwidth pinned by control_tds$Ht: a single candidate (see above).
        vkt <- opt_t[k]
      } else {
        if (opt_t[k] < max_dist_t) {
          up_cand_t <- Vt[Vt > opt_t[k]]
          up_t[k] <- if (length(up_cand_t)) tail(up_cand_t, 1) else opt_t[k]
        } else {
          up_t[k] <- max_dist_t
        }

        if (opt_t[k] > min_dist_t) {
          if (sum(Vt < opt_t[k]) > 0)
            down_t[k] <- head(Vt[Vt < opt_t[k]], 1)
          else
            down_t[k] <- min_dist_t
          ddown_t[k] <- max(min(down_t[down_t > min_dist_t], ddown_t[k], na.rm = TRUE),
                            min_dist_t, na.rm = TRUE)
        } else {
          ddown_t[k] <- down_t[k] <- min_dist_t
        }

        # --- Add extended candidates V5t ---
        if (!is.null(V5t)) {
          v5t_up <- V5t[V5t > up_t[k]]
          if (length(v5t_up) > 0) {
            up5t <- min(v5t_up)
            extra_vt <- V5t[V5t >= up5t]
            vkt <- sort(unique(c(max_dist_t, up_t[k], opt_t[k], down_t[k], ddown_t[k], extra_vt)), decreasing = TRUE)
          } else {
            vkt <- sort(unique(c(max_dist_t, up_t[k], opt_t[k], down_t[k], ddown_t[k])), decreasing = TRUE)
          }
        } else {
          vkt <- sort(unique(c(max_dist_t, up_t[k], opt_t[k], down_t[k], ddown_t[k])), decreasing = TRUE)
        }
      }

    } else {
      vkt <- NULL
    }

  })
}
