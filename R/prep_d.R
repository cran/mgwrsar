#' prep_d
#' to be documented
#' @usage prep_d(coords, NN, TP, extrapol = FALSE, ratio = 1, QP = NULL, Type = NULL)
#' @param coords A matrix with spatial coordinates
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP index of target points, default 1:n
#' @param extrapol  special mode for prediction unig GWR estimation, default FALSE
#' @param ratio numeric [0,1] ratio time/space for ordering indexG, default 1.
#' @param QP index of query points, default NULL
#' @param kernels kernels
#' @noRd
#' @return A list with precomputed matrices of distances and neighbours indexes
prep_d <- function(coords, NN, TP, Q = FALSE, extrapol = FALSE,
                   ratio = 1, QP = NULL, kernels = NULL, Type = NULL, stable_knn = FALSE) {

  if(length(TP) < nrow(coords)) stable_knn = FALSE
  dists <- list()
  if (is.null(QP)) QP <- 1:nrow(coords)
  TPc <- (1:nrow(coords))[-TP]   # currently unused, but kept

  # ---- Cyclic temporal KNN (as you implemented) ----
  knn_cyclic_time_modulo <- function(time_vec, k, query = NULL,
                                     period = 365, dedup = TRUE) {
    if (is.null(query)) query <- time_vec
    n <- length(time_vec)

    # 1) Augmented data (-P, 0, +P) + mapping to original indices
    time_aug <- c(time_vec - period, time_vec, time_vec + period)
    idx_map  <- rep.int(seq_len(n), 3L)

    # 2) KNN on the augmented 1D axis
    res <- nabor::knn(
      data  = matrix(time_aug,  ncol = 1),
      query = matrix(query,     ncol = 1),
      k = k * if (dedup) 3L else 1L
    )


    # 3) Remap returned indices to original space
    idx_orig <- matrix(idx_map[res$nn.idx],
                       nrow = nrow(res$nn.idx),
                       ncol = ncol(res$nn.idx))

    # 4) Cyclic distances
    tq   <- matrix(query,              nrow = nrow(idx_orig), ncol = ncol(idx_orig))
    tj   <- matrix(time_vec[idx_orig], nrow = nrow(idx_orig), ncol = ncol(idx_orig))
    dlin <- abs(tq - tj)
    dmod <- dlin %% period
    dcyc <- pmin(dmod, period - dmod)

    # 5) Deduplication by query
    if (dedup) {
      keep_idx <- matrix(NA_integer_, nrow = nrow(idx_orig), ncol = k)
      keep_dst <- matrix(NA_real_,    nrow = nrow(idx_orig), ncol = k)
      for (i in seq_len(nrow(idx_orig))) {
        ord <- order(dcyc[i, ], decreasing = FALSE)
        uniq_ord <- ord[!duplicated(idx_orig[i, ord])]
        sel <- head(uniq_ord, k)
        keep_idx[i, ] <- idx_orig[i, sel]
        keep_dst[i, ] <- dcyc[i, sel]
      }
      idx_orig <- keep_idx
      dcyc     <- keep_dst
    }

    list(nn.idx = idx_orig, nn.dists = dcyc)
  }

  ##================== 1) SPATIAL KNN ==================##
  if (Type %in% c("GD")) {
    if (extrapol) {
      nn <- nabor::knn(coords[-TP, 1:2],
                       query = coords[TP, 1:2],
                       k = min(c(NN, nrow(coords) - length(TP))))
      if (stable_knn) {
        knn_sorted <- .Call(
          "_mgwrsar_knn_stable_sort",
          nn$nn.dist,
          nn$nn.idx
        )
        nn$nn.dist <- knn_sorted$dist
        nn$nn.idx  <- knn_sorted$idx
      }
    } else {
      nn <- nabor::knn(coords[, 1:2],
                       query = coords[TP, 1:2],
                       k = NN)
      if (stable_knn) {
        knn_sorted <- .Call(
          "_mgwrsar_knn_stable_sort",
          nn$nn.dist,
          nn$nn.idx
        )
        nn$nn.dist <- knn_sorted$dist
        nn$nn.idx  <- knn_sorted$idx
      }
    }
    indexGS <- nn$nn.idx
    DS      <- nn$nn.dists
  }

  ##================== 2) TEMPORAL KNN =================##
  if (Type %in% c("T")) {
    if (Type == "T") mic <- 1 else mic <- 3

    if (extrapol) {
      nnt <- nabor::knn(coords[-TP, mic, drop = FALSE],
                        query = coords[TP, mic, drop = FALSE],
                        k = min(c(NN, nrow(coords) - length(TP))))
      if (stable_knn) {
        knn_sorted <- .Call(
          "_mgwrsar_knn_stable_sort",
          nnt$nn.dist,
          nnt$nn.idx
        )
        nnt$nn.dist <- knn_sorted$dist
        nnt$nn.idx  <- knn_sorted$idx
      }
    } else {
      nnt <- nabor::knn(coords[, mic, drop = FALSE],
                        query = coords[TP, mic, drop = FALSE],
                        k = NN)
      if (stable_knn) {
        knn_sorted <- .Call(
          "_mgwrsar_knn_stable_sort",
          nnt$nn.dist,
          nnt$nn.idx
        )
        nnt$nn.dist <- knn_sorted$dist
        nnt$nn.idx  <- knn_sorted$idx
      }
    }

    indexGT <- nnt$nn.idx
    DT      <- nnt$nn.dists

    # Retrieve potential cycling period from the last kernel
    cycling <- as.numeric(
      unlist(stringr::str_split(tail(kernels, 1), "_"))[3]
    )

    if (!is.na(cycling)) {
      # replace DT with modulo version
      if (extrapol) {
        nntm <- knn_cyclic_time_modulo(coords[-TP, mic],
                                       k       = min(c(NN, nrow(coords) - length(TP))),
                                       query   = coords[TP, mic],
                                       period = cycling)
        if (stable_knn) {
          knn_sorted <- .Call(
            "_mgwrsar_knn_stable_sort",
            nntm$nn.dist,
            nntm$nn.idx
          )
          nntm$nn.dist <- knn_sorted$dist
          nntm$nn.idx  <- knn_sorted$idx
        }

      } else {
        nntm <- knn_cyclic_time_modulo(coords[, mic],
                                       k       = NN,
                                       query   = coords[TP, mic],
                                       period = cycling)
        if (stable_knn) {
          knn_sorted <- .Call(
            "_mgwrsar_knn_stable_sort",
            nntm$nn.dist,
            nntm$nn.idx
          )
          nntm$nn.dist <- knn_sorted$dist
          nntm$nn.idx  <- knn_sorted$idx
        }
      }
      indexGT <- nntm$nn.idx
      DT      <- nntm$nn.dists  # DT = single time distance (modulo)
    }
  }

  ##================== 3) CONSTRUCTION dists / indexG =================##

  gdt_knn_scaled <- function(coords, Time, TP = NULL,
                             NN = 4000,
                             cycling = NULL,
                             extrapol = FALSE,
                             subsample = 4000,
                             seed = 123) {

    # ---- Fix: NA -> NULL for cycling ----
    if (length(cycling) == 1L && is.na(cycling)) cycling <- NULL

    stopifnot(nrow(coords) == length(Time))
    n <- nrow(coords)

    if (is.null(TP)) TP <- seq_len(n)

    ALL <- seq_len(n)
    TRAIN <- if (extrapol) setdiff(ALL, TP) else ALL

    coords_train <- coords[TRAIN, , drop = FALSE]
    time_train   <- Time[TRAIN]

    coords_query <- coords[TP, , drop = FALSE]
    time_query   <- Time[TP]

    # ---- True cyclic distance helper (modulo) ----
    cyclic_dist <- function(dt, cycling) {
      dt <- dt %% cycling
      pmin(dt, cycling - dt)
    }

    # =====================================================
    # 1. Estimate scaling factors (Option A)
    # =====================================================
    set.seed(seed)
    idx <- sample(TRAIN, min(subsample, length(TRAIN)))

    ## Spatial distances
    DS <- as.vector(dist(coords[idx, , drop = FALSE]))

    ## Temporal distances (correct modulo)
    DTmat <- abs(outer(Time[idx], Time[idx], "-"))
    if (!is.null(cycling)) DTmat <- cyclic_dist(DTmat, cycling)
    DT <- DTmat[upper.tri(DTmat)]

    scale_space <- sd(DS)
    scale_time  <- sd(DT)

    if (scale_space == 0 || scale_time == 0)
      stop("Zero variance in spatial or temporal distances.")

    # =====================================================
    # 2. Build scaled KNN space
    # =====================================================
    if (is.null(cycling)) {

      X_train <- cbind(
        coords_train[,1] / scale_space,
        coords_train[,2] / scale_space,
        time_train        / scale_time
      )

      X_query <- cbind(
        coords_query[,1] / scale_space,
        coords_query[,2] / scale_space,
        time_query        / scale_time
      )

    } else {

      theta_train <- 2 * pi * (time_train %% cycling) / cycling
      theta_query <- 2 * pi * (time_query %% cycling) / cycling

      X_train <- cbind(
        coords_train[,1] / scale_space,
        coords_train[,2] / scale_space,
        cos(theta_train) / scale_time,
        sin(theta_train) / scale_time
      )

      X_query <- cbind(
        coords_query[,1] / scale_space,
        coords_query[,2] / scale_space,
        cos(theta_query) / scale_time,
        sin(theta_query) / scale_time
      )
    }

    # =====================================================
    # 3. KNN search
    # =====================================================
    k_eff <- min(NN, nrow(X_train))

    kn <- nabor::knn(
      data  = X_train,
      query = X_query,
      k     = k_eff
    )

    indexG <- matrix(TRAIN[kn$nn.idx],
                     nrow = length(TP),
                     ncol = k_eff)

    # =====================================================
    # 4. Fast distance reconstruction (vectorised per row)
    # =====================================================
    DS_out <- matrix(NA_real_, length(TP), k_eff)
    DT_out <- matrix(NA_real_, length(TP), k_eff)

    # for (ii in seq_along(TP)) {
    #
    #   i  <- TP[ii]
    #   jj <- indexG[ii, ]
    #
    #   ## spatial
    #   dx <- coords[jj,1] - coords[i,1]
    #   dy <- coords[jj,2] - coords[i,2]
    #   DS_out[ii, ] <- sqrt(dx*dx + dy*dy)
    #
    #   ## temporal (correct modulo)
    #   dt_lin <- abs(Time[jj] - Time[i])
    #   DT_out[ii, ] <- if (is.null(cycling)) dt_lin else cyclic_dist(dt_lin, cycling)
    # }
    res <- .Call(
      "_mgwrsar_compute_DS_DT_cpp",
      coords,
      Time,
      indexG,
      !is.null(cycling),
      if (is.null(cycling)) 0 else cycling,
      PACKAGE = "mgwrsar"
    )

    list(
      indexG      = indexG,
      DS          = res$DS,
      DT          = res$DT
    )
  }

  if(Type == 'GDT'){
    cycling <- as.numeric(
      unlist(stringr::str_split(tail(kernels, 1), "_"))[3]
    )
    # if(NN==n+1){
    #   DS <- reord_M_R(DS, indexGS) # data order ; DS[i,j] indicates spatial distance between i and j
    #   DT = reord_M_R(DT, indexGT) # data order; DTM[i,j] indicates temporal modulo distance between i and j
    #   indexes <- get_index_mahalanobis_dual_apply(DS, DT)
    #   dists[['dist_t']] <- reord_M_R(DT, indexes$index_inverse)
    #   dists[['dist_s']] <- reord_M_R(DS, indexes$index_inverse)
    #   indexG <- indexes$index_order
    # } else{
      GDT<-gdt_knn_scaled(coords[,1:2],coords[,3],TP=TP,NN=NN, cycling = cycling,extrapol =extrapol)
      dists[['dist_t']] <- GDT$DT
      dists[['dist_s']] <- GDT$DS
      indexG <- GDT$indexG
    #}
  } else if(Type == 'GD'){
    indexG <- indexGS
    dists[['dist_s']] <- DS
  } else if(Type == 'T'){
    indexG <- indexGT
    dists[['dist_t']] <- DT
  }
  list(indexG = indexG, dists = dists)
}

