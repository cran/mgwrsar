update_bandwidth_candidates <- function(env = parent.frame()) {
  with(env, {
    # Safety: default values
    if (!exists("V5")) V5 <- NULL
    if (!exists("V5t")) V5t <- NULL

    # -----------------------------------------------------------------
    #  Spatial part: compute vks
    # -----------------------------------------------------------------
    if (opt[k] < max_dist)
      up[k] <- tail(V[V > opt[k]], 1)
    else
      up[k] <- max_dist

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

    # -----------------------------------------------------------------
    #  Temporal part: compute vkt (if opt_t and Vt exist)
    # -----------------------------------------------------------------
    if (exists("opt_t") && !is.null(opt_t) && exists("Vt") && !is.null(Vt)) {
      if (opt_t[k] < max_dist_t)
        up_t[k] <- tail(Vt[Vt > opt_t[k]], 1)
      else
        up_t[k] <- max_dist_t

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

    } else {
      vkt <- NULL
    }

  })
}
