#' prep_d
#' to be documented
#' @usage prep_d(coords, NN, TP, extrapol = FALSE, ratio = 1, QP = NULL, Type = NULL)
#' @param coords A matrix with spatial coordinates
#' @param NN Number of spatial Neighbours for kernels computations
#' @param TP index of target points, default 1:n
#' @param extrapol special mode for prediction unig GWR estimation, default FALSE
#' @param ratio numeric [0,1] ratio time/space for ordering indexG, default 1.
#' @param QP index of query points, default NULL
#' @param kernels kernels
#' @noRd
#' @return A list with precomputed matrices of distances and neighbours indexes

prep_d<-function (coords, NN, TP, Q = FALSE, extrapol = FALSE, ratio = 1,
                  QP = NULL, kernels = NULL, Type = NULL, stable_knn = FALSE)
{
  .mgwrsar_set_native_threads(1L)
  if (length(TP) < nrow(coords))
    stable_knn = FALSE
  dists <- list()
  if (is.null(QP))
    QP <- 1:nrow(coords)
  knn_cyclic_time_modulo <- function(time_vec, k, query = NULL,
                                     period = 365, dedup = TRUE) {
    if (is.null(query))
      query <- time_vec
    n <- length(time_vec)
    time_aug <- c(time_vec - period, time_vec, time_vec +
                    period)
    idx_map <- rep.int(seq_len(n), 3L)
    res <- nabor::knn(data = matrix(time_aug, ncol = 1),
                      query = matrix(query, ncol = 1), k = k * if (dedup)
                        3L
                      else 1L)
    idx_orig <- matrix(idx_map[res$nn.idx], nrow = nrow(res$nn.idx),
                       ncol = ncol(res$nn.idx))
    tq <- matrix(query, nrow = nrow(idx_orig), ncol = ncol(idx_orig))
    tj <- matrix(time_vec[idx_orig], nrow = nrow(idx_orig),
                 ncol = ncol(idx_orig))
    dlin <- abs(tq - tj)
    dmod <- dlin%%period
    dcyc <- pmin(dmod, period - dmod)
    if (dedup) {
      keep_idx <- matrix(NA_integer_, nrow = nrow(idx_orig),
                         ncol = k)
      keep_dst <- matrix(NA_real_, nrow = nrow(idx_orig),
                         ncol = k)
      for (i in seq_len(nrow(idx_orig))) {
        ord <- order(dcyc[i, ], decreasing = FALSE)
        uniq_ord <- ord[!duplicated(idx_orig[i, ord])]
        sel <- head(uniq_ord, k)
        keep_idx[i, ] <- idx_orig[i, sel]
        keep_dst[i, ] <- dcyc[i, sel]
      }
      idx_orig <- keep_idx
      dcyc <- keep_dst
    }
    list(nn.idx = idx_orig, nn.dists = dcyc)
  }
  if (Type %in% c("GD")) {
    if (extrapol) {
      nn <- nabor::knn(coords[-TP, 1:2], query = coords[TP,
                                                        1:2], k = min(c(NN, nrow(coords) - length(TP))))
    }
    else {
      nn <- nabor::knn(coords[, 1:2], query = coords[TP,
                                                     1:2], k = NN)
    }
    if (stable_knn) {
      knn_sorted <- .Call("_mgwrsar_knn_stable_sort", nn$nn.dist,
                          nn$nn.idx)
      nn$nn.dist <- knn_sorted$dist
      nn$nn.idx <- knn_sorted$idx
    }
    indexGS <- nn$nn.idx
    DS <- nn$nn.dists
  }
  if (Type %in% c("T")) {
    mic <- if (Type == "T")
      1
    else 3
    col_idx <- if (ncol(coords) >= 3)
      3
    else 1
    if (Type == "T")
      col_idx <- 1
    if (extrapol) {
      nnt <- nabor::knn(coords[-TP, col_idx, drop = FALSE],
                        query = coords[TP, col_idx, drop = FALSE], k = min(c(NN,
                                                                             nrow(coords) - length(TP))))
    }
    else {
      nnt <- nabor::knn(coords[, col_idx, drop = FALSE],
                        query = coords[TP, col_idx, drop = FALSE], k = NN)
    }
    if (stable_knn) {
      knn_sorted <- .Call("_mgwrsar_knn_stable_sort", nnt$nn.dist,
                          nnt$nn.idx)
      nnt$nn.dist <- knn_sorted$dist
      nnt$nn.idx <- knn_sorted$idx
    }
    indexGT <- nnt$nn.idx
    DT <- nnt$nn.dists
    cycling <- as.numeric(unlist(stringr::str_split(tail(kernels,
                                                         1), "_"))[3])
    if (!is.na(cycling)) {
      if (extrapol) {
        nntm <- knn_cyclic_time_modulo(coords[-TP, col_idx],
                                       k = min(c(NN, nrow(coords) - length(TP))),
                                       query = coords[TP, col_idx], period = cycling)
      }
      else {
        nntm <- knn_cyclic_time_modulo(coords[, col_idx],
                                       k = NN, query = coords[TP, col_idx], period = cycling)
      }
      if (stable_knn) {
        knn_sorted <- .Call("_mgwrsar_knn_stable_sort",
                            nntm$nn.dists, nntm$nn.idx)
        nntm$nn.dists <- knn_sorted$dist
        nntm$nn.idx <- knn_sorted$idx
      }
      indexGT <- nntm$nn.idx
      DT <- nntm$nn.dists
    }
  }
  gdt_knn_scaled <- function(coords, Time, TP = NULL, NN = 4000,
                             cycling = NULL, extrapol = FALSE, subsample = 4000, seed = 123) {
    if (length(cycling) == 1L && is.na(cycling))
      cycling <- NULL
    stopifnot(nrow(coords) == length(Time))
    n <- nrow(coords)
    if (is.null(TP))
      TP <- seq_len(n)
    ALL <- seq_len(n)
    TRAIN <- if (extrapol)
      setdiff(ALL, TP)
    else ALL
    coords_train <- coords[TRAIN, , drop = FALSE]
    time_train <- Time[TRAIN]
    coords_query <- coords[TP, , drop = FALSE]
    time_query <- Time[TP]
    cyclic_dist <- function(dt, cycling) {
      dt <- dt%%cycling
      pmin(dt, cycling - dt)
    }
    rng_state <- .mgwrsar_rng_save(); on.exit(.mgwrsar_rng_restore(rng_state), add = TRUE)
    set.seed(seed)
    idx <- sample(TRAIN, min(subsample, length(TRAIN)))
    DS_samp <- as.vector(dist(coords[idx, , drop = FALSE]))
    DTmat <- abs(outer(Time[idx], Time[idx], "-"))
    if (!is.null(cycling))
      DTmat <- cyclic_dist(DTmat, cycling)
    DT_samp <- DTmat[upper.tri(DTmat)]
    scale_space <- sd(DS_samp)
    scale_time <- sd(DT_samp)
    if (scale_space == 0)
      scale_space <- 1
    if (scale_time == 0)
      scale_time <- 1
    if (is.null(cycling)) {
      X_train <- cbind(coords_train[, 1]/scale_space, coords_train[,
                                                                   2]/scale_space, time_train/scale_time)
      X_query <- cbind(coords_query[, 1]/scale_space, coords_query[,
                                                                   2]/scale_space, time_query/scale_time)
    }
    else {
      theta_train <- 2 * pi * (time_train%%cycling)/cycling
      theta_query <- 2 * pi * (time_query%%cycling)/cycling
      X_train <- cbind(coords_train[, 1]/scale_space, coords_train[,
                                                                   2]/scale_space, cos(theta_train)/scale_time,
                       sin(theta_train)/scale_time)
      X_query <- cbind(coords_query[, 1]/scale_space, coords_query[,
                                                                   2]/scale_space, cos(theta_query)/scale_time,
                       sin(theta_query)/scale_time)
    }
    k_eff <- min(NN, nrow(X_train))
    kn <- nabor::knn(data = X_train, query = X_query, k = k_eff)
    indexG <- matrix(TRAIN[kn$nn.idx], nrow = length(TP),
                     ncol = k_eff)
    n_tp <- length(TP)
    X_neigh <- matrix(coords[indexG, 1], nrow = n_tp, ncol = k_eff)
    Y_neigh <- matrix(coords[indexG, 2], nrow = n_tp, ncol = k_eff)
    X_targ <- coords[TP, 1]
    Y_targ <- coords[TP, 2]
    dX <- X_neigh - X_targ
    dY <- Y_neigh - Y_targ
    DS_out <- sqrt(dX^2 + dY^2)
    T_neigh <- matrix(Time[indexG], nrow = n_tp, ncol = k_eff)
    T_targ <- Time[TP]
    dT_raw <- abs(T_neigh - T_targ)
    if (!is.null(cycling)) {
      dT_mod <- dT_raw%%cycling
      DT_out <- pmin(dT_mod, cycling - dT_mod)
    }
    else {
      DT_out <- dT_raw
    }
    list(indexG = indexG, DS = DS_out, DT = DT_out)
  }
  if (Type == "GDT") {
    cycling <- as.numeric(unlist(stringr::str_split(tail(kernels,
                                                         1), "_"))[3])
    GDT <- gdt_knn_scaled(coords[, 1:2], coords[, 3], TP = TP,
                          NN = NN, cycling = cycling, extrapol = extrapol)
    dists[["dist_t"]] <- GDT$DT
    dists[["dist_s"]] <- GDT$DS
    indexG <- GDT$indexG
  }
  else if (Type == "GD") {
    indexG <- indexGS
    dists[["dist_s"]] <- DS
  }
  else if (Type == "T") {
    indexG <- indexGT
    dists[["dist_t"]] <- DT
  }
  list(indexG = indexG, dists = dists)
}
