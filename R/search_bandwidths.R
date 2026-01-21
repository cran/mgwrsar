#' Bandwidth Selection via Multi-round Grid Search based on AICc
#'
#' @description
#' This function selects the optimal bandwidths (spatial and/or temporal) for GWR or MGWR models
#' by minimizing the Corrected Akaike Information Criterion (AICc). It uses a multi-round
#' grid search approach (coarse-to-fine) to efficiently narrow down the optimal parameter space
#' before optionally applying a golden section search for final refinement.
#'
#' @param formula A formula object specifying the model (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing the variables in the model.
#' @param coords A matrix or data frame of coordinates (2 columns for spatial, 1 for temporal, or more for GDT).
#' @param fixed_vars A character vector indicating the names of variables with spatially stationary (fixed) coefficients.
#' Default is \code{NULL} (all coefficients are varying).
#' @param kernels A character vector specifying the kernel types for spatial and temporal components
#' (e.g., \code{c("gauss", "gauss")}).
#' @param Model A character string specifying the model type. Options include "GWR", "MGWR", "OLS", "SAR", etc.
#' Default is "GWR".
#' @param control A named list of extra control arguments passed to the \code{MGWRSAR} function
#' (e.g., \code{adaptive}, \code{NN}, \code{Z}).
#' @param hs_range A numeric vector of length 2 defining the lower and upper bounds for the spatial bandwidth search.
#' @param ht_range A numeric vector of length 2 defining the lower and upper bounds for the temporal bandwidth search.
#' Set to \code{NULL} for spatial-only models.
#' @param n_seq An integer specifying the number of bandwidth candidates to test per dimension in each round.
#' @param ncore An integer specifying the number of CPU cores to use for parallel processing.
#' Default is \code{parallel::detectCores() - 1}.
#' @param n_rounds An integer specifying the number of grid search rounds (zooming steps). Default is 3.
#' @param refine Logical. If \code{TRUE}, a final optimization step using golden section search is performed
#' around the best candidate found. Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages to the console.
#' @param show_progress Logical. If \code{TRUE}, displays a progress bar during computation.
#' @param tol A numeric vector of length 2 (or 1) specifying the
#'  tolerance for spatial and temporal bandwidths. If \code{NULL}, it
#'  is calculated automatically based on the range.
#' @param parallel_method Parallelization method ("auto", "fork", "socket"); "auto" selects the best available backend depending on
#' the OS, other values run sequentially.
#'
#' @details
#' The function performs a grid search over \code{n_rounds}. In the first round, it tests \code{n_seq}
#' candidates linearly or geometrically spaced within the provided ranges. In subsequent rounds,
#' the search range is narrowed around the best candidate from the previous round.
#'
#'
#' @return A list containing:
#' \item{mode}{The mode of optimization ("spatial_only" or "spatio-temporal").}
#' \item{results}{A list of data frames containing the results (bandwidths and AICc) for each round.}
#' \item{best}{A data frame row corresponding to the best parameter combination found.}
#' \item{refined}{The result of the refinement step (if \code{refine = TRUE}).}
#' \item{best_model}{The final \code{mgwrsar} model object fitted with the optimal bandwidths.}
#' \item{ctime}{The total computation time.}
#'
#' @seealso \code{\link{MGWRSAR}}
#' @export
search_bandwidths<- function(
    formula, data, coords, fixed_vars = NULL,
    kernels = c("gauss","gauss"),
    Model = "GWR",
    control = list(),
    hs_range = NULL,
    ht_range = NULL,
    n_seq = 10,
    ncore = 1,
    n_rounds =0,
    refine = TRUE,
    verbose = FALSE,
    show_progress = FALSE,
    tol = NULL,
    parallel_method = "auto"
){

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  # -------------------------------------------------------------------------
  # 0) Sanity checks on required ranges (by Type)
  # -------------------------------------------------------------------------
  if (is.null(control$Type)) control$Type <- "GD"
  if (is.null(control$TP)) control$TP<-1:nrow(data)

  if (control$Type %in% "GDT" && is.null(hs_range) && is.null(ht_range)) {
    stop(
      "For spatio-temporal models (Type = 'GDT'), both `hs_range` (spatial) and ",
      "`ht_range` (temporal) must be provided."
    )
  }
  if (control$Type %in% c("GD", "GDT") && is.null(hs_range)) {
    stop("Missing `hs_range`: spatial bandwidth range is required for Type = '", control$Type, "'.")
  }
  if (control$Type %in% c("T", "GDT") && is.null(ht_range)) {
    stop("Missing `ht_range`: temporal bandwidth range is required for Type = '", control$Type, "'.")
  }

  # -------------------------------------------------------------------------
  # 1) Normalize parallel control for MGWRSAR (your helper)
  # -------------------------------------------------------------------------
  control <- .mgwrsar_normalize_parallel_control(control, context = "MGWRSAR")

  # -------------------------------------------------------------------------
  # 2) Criterion utilities (shared by grid + refine decision)
  # -------------------------------------------------------------------------
  .get_criterion_value <- function(mod, crit) {
    if (is.null(crit) || !nzchar(crit)) crit <- "AICc"
    slots <- slotNames(mod)

    if (crit %in% c("CV","CVtp")) {
      if (crit == "CV"   && "RMSE"   %in% slots) return(as.numeric(mod@RMSE))
      if (crit == "CVtp" && "RMSEtp" %in% slots) return(as.numeric(mod@RMSEtp))
      return(NA_real_)
    }

    if (crit %in% slots) return(as.numeric(slot(mod, crit)))
    if ("AICc" %in% slots) return(as.numeric(mod@AICc))
    if ("AIC"  %in% slots) return(as.numeric(mod@AIC))
    NA_real_
  }

  # -------------------------------------------------------------------------
  # 3) Internal golden search (legacy code, but kept internal)
  #    -> IMPORTANT: safe_eval() already uses control$criterion
  # -------------------------------------------------------------------------

  golden_search_bandwidth_in <- function(
    formula, Ht = NULL, data, coords, fixed_vars = NULL, kernels,
    Model = "GWR", control, lower.bound, upper.bound,
    tolerance = 1e-6,
    ncore = 1,
    show_progress = TRUE
  ) {

    ptm <- proc.time()
    `%||%` <- function(x, y) if (!is.null(x)) x else y

    get_SFdata()

    mf <- model.frame(formula, data)
    data <- data[, names(mf), drop = FALSE]

    acc <- new.env(parent = emptyenv())
    acc$isolated_hits <- 0L

    # -- control setup
    control_o <- control
    control$SE     <- FALSE
    control$get_s  <- FALSE
    control$get_Rk <- FALSE
    if (is.null(control$Type)) control$Type <- "GD"

    # -- jitter coords
    if (control$Type %in% c("GD","GDT")) coords <- make_unique_by_structure(coords)
    if (control$Type %in% c("GDT","T"))  control$Z <- make_unique_by_structure(control$Z)

    adaptive <- isTRUE(control$adaptive[1])
    golden.ratio <- 2 / (sqrt(5) + 1)

    # -- starting points
    x1 <- upper.bound - golden.ratio * (upper.bound - lower.bound)
    x2 <- lower.bound + golden.ratio * (upper.bound - lower.bound)
    if (adaptive) { x1 <- floor(x1); x2 <- ceiling(x2) }

    # -- distances computation (heavy)
    if (is.null(control$dists)) {
      if (is.null(control$NN)) control$NN <- nrow(data)
      if (is.null(control$TP)) control$TP <- 1:nrow(data)

      S <- if (control$Type == "T") as.matrix(control$Z, ncol = 1)
      else if (control$Type == "GDT") as.matrix(cbind(coords, control$Z), ncol = 3)
      else as.matrix(coords)

      stage1 <- prep_d(coords = S, NN = control$NN, TP = control$TP,
                       kernels = kernels, Type = control$Type)
      control$indexG <- stage1$indexG
      control$dists  <- stage1$dists
    }
    control_o <- control

    # --- safe_eval: uses control$criterion coherently ---
    safe_eval <- function(H) {
      crit <- control$criterion %||% "AICc"

      if (crit %in% c("CV","CVtp")) {
        control2 <- control
        control2$isgcv <- TRUE

        mod <- tryCatch(
          MGWRSAR(formula = formula, data = data, coords = coords, fixed_vars = fixed_vars,
                  kernels = kernels, H = H, Model = Model, control = control2),
          error = function(e) NULL
        )
        if (is.null(mod)) return(Inf)

        slots <- slotNames(mod)
        if (crit == "CV"   && "RMSE"   %in% slots) return(as.numeric(mod@RMSE))
        if (crit == "CVtp" && "RMSEtp" %in% slots) return(as.numeric(mod@RMSEtp))
        return(Inf)

      } else {
        control_AIC <- control
        control_AIC$get_ts <- TRUE

        mod <- tryCatch(
          MGWRSAR(formula = formula, data = data, coords = coords, fixed_vars = fixed_vars,
                  kernels = kernels, H = H, Model = Model, control = control_AIC),
          error = function(e) NULL
        )
        if (is.null(mod)) return(Inf)

        val <- .get_criterion_value(mod, crit)
        if (!is.finite(val)) return(Inf)
        return(val)
      }
    }

    # --- parallel inside golden (kept, but safe) ---
    cl <- NULL
    use_parallel <- (ncore > 1) &&
      requireNamespace("doParallel", quietly = TRUE) &&
      requireNamespace("foreach", quietly = TRUE)

    if (use_parallel) {
      sysname <- Sys.info()[["sysname"]]
      if (sysname == "Linux") {
        doParallel::registerDoParallel(cores = ncore)
      } else {
        cl <- parallel::makeCluster(ncore)
        doParallel::registerDoParallel(cl)
        parallel::clusterEvalQ(cl, {
          suppressPackageStartupMessages(library(mgwrsar))
          NULL
        })
        parallel::clusterExport(
          cl,
          varlist = c("formula","data","coords","fixed_vars","kernels","Model","control","safe_eval"),
          envir = environment()
        )
      }
      on.exit({
        if (!is.null(cl)) parallel::stopCluster(cl)
        doParallel::stopImplicitCluster()
      }, add = TRUE)
    }

    eval_pair_parallel <- function(h1, h2) {
      H_list <- list(c(h1, Ht), c(h2, Ht))
      if (use_parallel) {
        foreach::foreach(h = H_list, .combine = c, .packages = c("mgwrsar")) %dopar% safe_eval(h)
      } else {
        sapply(H_list, safe_eval)
      }
    }

    f12 <- eval_pair_parallel(x1, x2)
    f1 <- f12[1]; f2 <- f12[2]

    iteration <- 0L
    if (adaptive) tolerance <- 1
    x1_p <- 1; x2_p <- 2
    max_iter <- 50L

    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }

    while (abs(upper.bound - lower.bound) > tolerance &&
           ((adaptive && abs(x1 - x2) > 1 && (x1_p != x1 || x2_p != x2)) || !adaptive)) {

      iteration <- iteration + 1L
      if (show_progress) setTxtProgressBar(pb, iteration)
      x1_p <- x1; x2_p <- x2

      if (f2 > f1) {
        upper.bound <- x2
        x2 <- x1; f2 <- f1
        x1 <- upper.bound - golden.ratio * (upper.bound - lower.bound)
        if (adaptive) x1 <- floor(x1)
        f12 <- eval_pair_parallel(x1, x2)
        f1 <- f12[1]; f2 <- f12[2]
      } else {
        lower.bound <- x1
        x1 <- x2; f1 <- f2
        x2 <- lower.bound + golden.ratio * (upper.bound - lower.bound)
        if (adaptive) x2 <- ceiling(x2)
        f12 <- eval_pair_parallel(x1, x2)
        f1 <- f12[1]; f2 <- f12[2]
      }

      if (iteration >= max_iter) break
    }

    res <- (lower.bound + upper.bound) / 2
    if (adaptive) {
      candidates <- unique(c(x1, x2, lower.bound, upper.bound))
      scores <- vapply(candidates, function(h) safe_eval(c(h, Ht)), numeric(1))
      res <- candidates[which.min(scores)]
      objective <- min(scores)
    } else {
      objective <- safe_eval(c(res, Ht))
    }

    final_model <- MGWRSAR(
      formula, data, coords, fixed_vars, kernels,
      H = c(res, Ht), Model = Model, control = control_o
    )

    ctime <- (proc.time() - ptm)[3]

    list(minimum = res, objective = objective, model = final_model, ctime = ctime)
  }

  # -------------------------------------------------------------------------
  # 4) Internal checks + data prep
  # -------------------------------------------------------------------------
  if (exists("get_SFdata")) get_SFdata()
  acc <- new.env(parent = emptyenv())
  acc$isolated_hits <- 0L

  n <- nrow(coords)
  set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")

  if (control$Type %in% c("GD","GDT")) coords <- make_unique_by_structure(coords)
  if (control$Type %in% c("GDT","T"))  control$Z <- make_unique_by_structure(control$Z)

  if (is.null(control$NN)) control$NN <- n
  mf <- model.frame(formula, data)
  data <- data[, names(mf), drop = FALSE]

  # -------------------------------------------------------------------------
  # 5) Heavy distance structures (unchanged)
  # -------------------------------------------------------------------------
  HEAVY <- list(
    data = data, coords = coords, fixed_vars = fixed_vars, Z = control$Z,
    main = list(indexG = NULL, dists = NULL),
    VT   = list(indexG = NULL, dists = NULL),
    VD   = list(indexG = NULL, dists = NULL)
  )

  VTcontrol <- list(); VDcontrol <- list()

  if (control$Type == "GDT") {
    G <- prep_d(coords = as.matrix(cbind(coords, control$Z)), NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
    HEAVY$main$indexG <- G$indexG; HEAVY$main$dists <- G$dists

    G_vt <- prep_d(coords = as.matrix(control$Z, ncol = 1), NN = control$NN, TP = 1:n, kernels = kernels[2], Type = "T")
    HEAVY$VT$indexG <- G_vt$indexG; HEAVY$VT$dists <- G_vt$dists
    VTcontrol <- control; VTcontrol$Type <- "T"; VTcontrol$adaptive <- FALSE
    VTcontrol$indexG <- NULL; VTcontrol$dists <- NULL; VTcontrol$Z <- NULL

    G_vd <- prep_d(coords = coords, NN = control$NN, TP =control$TP, kernels = kernels[1], Type = "GD")
    HEAVY$VD$indexG <- G_vd$indexG; HEAVY$VD$dists <- G_vd$dists
    VDcontrol <- control; VDcontrol$Type <- "GD"; VDcontrol$adaptive <- control$adaptive[1]
    VDcontrol$indexG <- NULL; VDcontrol$dists <- NULL; VDcontrol$Z <- NULL

    control <- unserialize(serialize(control, NULL))
    VTcontrol <- unserialize(serialize(VTcontrol, NULL))
    VDcontrol <- unserialize(serialize(VDcontrol, NULL))
  } else if (control$Type == "GD") {
    G <- prep_d(coords = coords, NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
    HEAVY$main$indexG <- G$indexG; HEAVY$main$dists <- G$dists
  }

  control$indexG <- NULL; control$dists <- NULL; control$Z <- NULL

  # -------------------------------------------------------------------------
  # 6) Parallel setup for GRID only (unchanged)
  # -------------------------------------------------------------------------
  if (ncore < 1) ncore <- 1
  sysname <- Sys.info()[["sysname"]]

  Sys.setenv(
    OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
  }

  use_parallel <- FALSE
  cl <- NULL

  if (ncore > 1 && requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
    use_parallel <- TRUE

    if (parallel_method == "auto") {
      type <- if (sysname != "Windows") "FORK" else "PSOCK"
    } else if (parallel_method == "fork") {
      if (sysname == "Windows") { warning("Fork not supported on Windows, switching to socket"); type <- "PSOCK" }
      else type <- "FORK"
    } else if (parallel_method == "socket") {
      type <- "PSOCK"
    } else {
      type <- "SEQ"
      use_parallel <- FALSE
    }

    if (use_parallel) {
      if (type == "FORK") {
        doParallel::registerDoParallel(cores = ncore)
        if (verbose) message(sprintf("Parallel backend: doParallel (Forking) with %d cores", ncore))
      } else {
        cl <- parallel::makeCluster(ncore)
        doParallel::registerDoParallel(cl)
        parallel::clusterEvalQ(cl, library(mgwrsar))
        if (verbose) message(sprintf("Parallel backend: doParallel (Socket) with %d cores", ncore))
      }
    }
  } else {
    foreach::registerDoSEQ()
  }

  on.exit({
    if (!is.null(cl)) parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    closeAllConnections()
  }, add = TRUE)

  # -------------------------------------------------------------------------
  # 7) Preparation
  # -------------------------------------------------------------------------
  start_time <- Sys.time()
  control_o <- control
  control$SE <- FALSE; control$get_s <- FALSE
  adaptive <- if ("adaptive" %in% names(control)) control$adaptive else c(FALSE, FALSE)
  if (length(tol) == 1) tol <- rep(tol, 2)

  round_timer <- function(expr, label = "", n_tests = 0) {
    t0 <- Sys.time()
    if (show_progress) message(sprintf("... [Round %s] Computing %d models ...", label, n_tests))
    val <- eval(expr)
    t1 <- Sys.time()
    dt <- as.numeric(difftime(t1, t0, units = "mins"))
    if (show_progress || verbose) message(sprintf("    [Round %s] Done in %.2f min.", label, dt))
    val
  }

  make_seq <- function(r) {
    if (r[1] <= 1e-9) seq(max(0, r[1]), r[2], length.out = n_seq)
    else exp(seq(log(r[1]), log(r[2]), length.out = n_seq))
  }
  get_subrange <- function(s, v, w = 2) {
    idx <- which.min(abs(s - v))
    c(s[max(1, floor(idx - w))], s[min(length(s), ceiling(idx + w))])
  }
  mid_subrange <- function(s, v) {
    r2 <- get_subrange(s, v, 2); r1 <- get_subrange(s, v, 1)
    c(mean(c(r2[1], r1[1])), mean(c(r2[2], r1[2])))
  }
  round_adaptive <- function(v, a, t, n) {
    v <- if (t > 0) round(v / t) * t else v
    if (isTRUE(a)) v <- v[v >= 1] else v <- v[v >= 0]
    v <- sort(unique(v))
    if (length(v) < min(3, n)) return(v)
    v
  }
  auto_tol <- function(r) {
    w <- max(r, na.rm = TRUE) - min(r, na.rm = TRUE)
    if (w <= 0) return(0)
    raw <- max(0.001 * w, 0.025 * w)
    k <- floor(log10(raw))
    m <- raw / (10^k)
    (if (m < 1.5) 1 else if (m < 3.5) 2 else if (m < 7.5) 5 else 10) * 10^k
  }

  if (missing(tol) || is.null(tol)) {
    tol_s <- if (!is.null(control$adaptive) && isTRUE(control$adaptive[1])) 1 else auto_tol(hs_range)
    tol_t <- if (length(ht_range) > 1) auto_tol(ht_range) else tol_s
    tol <- c(tol_s, tol_t)
  }

  hs_seq0 <- round_adaptive(make_seq(hs_range), adaptive[1], tol[1], n_seq)
  ht_seq0 <- if (!is.null(ht_range)) round_adaptive(make_seq(ht_range), adaptive[2], tol[2], n_seq) else NULL

  # -------------------------------------------------------------------------
  # 8) Worker: criterion-aware scoring
  # -------------------------------------------------------------------------
  compute_single_H <- function(i, H_grid, has_time, max_dist, max_dist_t,
                               formula, kernels, control, VTcontrol, VDcontrol, HEAVY) {

    if (has_time) { h_s <- H_grid[i, 1]; h_t <- H_grid[i, 2] }
    else { h_s <- H_grid[i, "h_s"]; h_t <- NA_real_ }

    control$indexG <- HEAVY$main$indexG; control$dists <- HEAVY$main$dists; control$Z <- HEAVY$Z
    if (has_time) {
      VTcontrol$indexG <- HEAVY$VT$indexG; VTcontrol$dists <- HEAVY$VT$dists; VTcontrol$Z <- HEAVY$Z
      VDcontrol$indexG <- HEAVY$VD$indexG; VDcontrol$dists <- HEAVY$VD$dists; VDcontrol$Z <- HEAVY$Z
    }

    mod <- tryCatch({
      if (!has_time) {
        if (h_s >= max_dist) {
          MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = NULL, Model = "OLS", control = control)
        } else {
          MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = h_s, Model = "GWR", control = control)
        }
      } else {
        if (h_s < max_dist && h_t < max_dist_t) {
          MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = c(h_s, h_t), Model = "GWR", control = control)
        } else if (h_s >= max_dist && h_t >= max_dist_t) {
          MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = NULL, Model = "OLS", control = control)
        } else if (h_s < max_dist && h_t >= max_dist_t) {
          MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels[1], H = c(h_s, NULL), Model = "GWR", control = VDcontrol)
        } else {
          MGWRSAR(formula, HEAVY$data, coords = as.matrix(HEAVY$Z, ncol = 1), fixed_vars = HEAVY$fixed_vars, kernels = kernels[2], H = h_t, Model = "GWR", control = VTcontrol)
        }
      }
    }, error = function(e) NULL)

    if (is.null(mod)) return(c(h_s = h_s, h_t = h_t, score = NA_real_, isolated = 0))

    crit <- control$criterion %||% "AICc"
    score <- .get_criterion_value(mod, crit)

    iso <- if ("isolated_idx" %in% slotNames(mod) && length(mod@isolated_idx) > 0) 1L else 0L
    c(h_s = h_s, h_t = h_t, score = score, isolated = iso)
  }

  # -------------------------------------------------------------------------
  # 9) Grid rounds: always minimize "score"
  # -------------------------------------------------------------------------
  run_rounds <- function(is_spatial_only) {
    hs_seq <- hs_seq0
    ht_seq <- if (!is_spatial_only) ht_seq0 else NULL
    results <- list(); best <- NULL

    if (n_rounds == 0) return(NULL)

    for (round_id in seq_len(min(n_rounds, 3))) {

      if (is_spatial_only) {
        if (round_id > 1) {
          hs_sub <- if (round_id == 2) mid_subrange(hs_seq, best$h_s) else get_subrange(hs_seq, best$h_s, 1)
          hs_seq <- round_adaptive(make_seq(hs_sub), adaptive[1], tol[1], n_seq)
        }
        H_grid <- data.frame(h_s = hs_seq)
        max_dist <- max(H_grid$h_s, na.rm = TRUE); max_dist_t <- NULL
      } else {
        if (round_id > 1) {
          hs_sub <- if (round_id == 2) mid_subrange(hs_seq, best$h_s) else get_subrange(hs_seq, best$h_s, 1)
          ht_sub <- if (round_id == 2) mid_subrange(ht_seq, best$h_t) else get_subrange(ht_seq, best$h_t, 1)
          hs_seq <- round_adaptive(make_seq(hs_sub), adaptive[1], tol[1], n_seq)
          ht_seq <- round_adaptive(make_seq(ht_sub), adaptive[2], tol[2], n_seq)
        }
        H_grid <- expand.grid(h_s = hs_seq, h_t = ht_seq)
        max_dist <- max(H_grid$h_s, na.rm = TRUE); max_dist_t <- max(H_grid$h_t, na.rm = TRUE)
      }

      if (verbose) {
        message(sprintf(">>> Round %d <<<", round_id))
        message(paste("h_s_seq =", paste(sprintf("%.6g", hs_seq), collapse = ", ")))
        if (!is_spatial_only) message(paste("h_t_seq =", paste(sprintf("%.6g", ht_seq), collapse = ", ")))
      }

      res_round <- round_timer({

        op <- if (use_parallel) foreach::`%dopar%` else foreach::`%do%`

        iter_obj <- foreach::foreach(
          i = 1:nrow(H_grid),
          .combine = rbind,
          .inorder = FALSE,
          .packages = c("mgwrsar")
        )

        out_mat <- op(iter_obj, {
          compute_single_H(
            i, H_grid, has_time = !is_spatial_only,
            max_dist = max_dist, max_dist_t = max_dist_t,
            formula = formula, kernels = kernels,
            control = control, VTcontrol = VTcontrol, VDcontrol = VDcontrol,
            HEAVY = HEAVY
          )
        })

        if (is.vector(out_mat) && !is.matrix(out_mat)) {
          out_mat <- matrix(out_mat, nrow = 1)
          colnames(out_mat) <- c("h_s", "h_t", "score", "isolated")
        }
        as.data.frame(out_mat)

      }, label = round_id, n_tests = nrow(H_grid))

      res_round$score <- as.numeric(res_round$score)
      if ("isolated" %in% names(res_round)) acc$isolated_hits <- acc$isolated_hits + sum(res_round$isolated, na.rm = TRUE)

      valid_idx <- which(is.finite(res_round$score))
      if (length(valid_idx) > 0) {
        best <- res_round[valid_idx[which.min(res_round$score[valid_idx])], , drop = FALSE]
      } else {
        best <- res_round[1, , drop = FALSE]
      }

      results[[round_id]] <- res_round

      if (verbose) {
        crit_lab <- control$criterion %||% "AICc"
        msg <- sprintf("Best h_s=%.6g", best$h_s)
        if (!is_spatial_only) msg <- paste0(msg, sprintf(", h_t=%.6g", best$h_t))
        message(paste0("   -> ", msg, sprintf(" (%s=%.6g)\n", crit_lab, best$score)))
      }
    }

    list(results = results, best = best, hs_seq = hs_seq, ht_seq = ht_seq)
  }

  # -------------------------------------------------------------------------
  # 10) EXECUTION
  # -------------------------------------------------------------------------
  is_spatial_only <- is.null(ht_range)
  control_o$indexG <- HEAVY$main$indexG; control_o$dists <- HEAVY$main$dists; control_o$Z <- HEAVY$Z

  if (n_rounds == 0) {
    do_p <- use_parallel
    #message("[Running] Direct Golden Search...")

    if (is_spatial_only) {
      refined <- tryCatch(
        golden_search_bandwidth_in(
          formula = formula, Ht = NULL,
          data = HEAVY$data, coords = HEAVY$coords, fixed_vars = HEAVY$fixed_vars,
          kernels = kernels, Model = Model,
          ncore = ifelse(do_p, ncore, 1),
          show_progress = show_progress,
          control = control_o,
          lower.bound = min(hs_seq0), upper.bound = max(hs_seq0),
          tolerance = tol[1]
        ),
        error = function(e) NULL
      )
    } else {
      refined <- tryCatch(
        golden_search_2d_bandwidth(
          formula = formula, data = HEAVY$data, coords = HEAVY$coords, fixed_vars = HEAVY$fixed_vars,
          kernels = kernels, Model = Model, control = control_o,
          lower.bound.space = min(hs_seq0), upper.bound.space = max(hs_seq0),
          lower.bound.time  = min(ht_seq0), upper.bound.time  = max(ht_seq0),
          tolerance_s = tol[1], tolerance_t = tol[2]
        ),
        error = function(e) NULL
      )
    }

    best_model <- if (!is.null(refined)) refined$model else NULL
    return(list(
      mode = if (is_spatial_only) "spatial" else "spatio-temporal",
      criterion = control$criterion %||% "AICc",
      results = NULL,
      best = NULL,
      refined = refined,
      objective=refined$objective,
      minimum=refined$minimum,
      best_model = best_model,
      ctime = refined$ctime
    ))
  }

  run_res <- run_rounds(is_spatial_only)
  best <- run_res$best

  # -------------------------------------------------------------------------
  # 11) Optional refinement (golden) using SAME criterion
  # -------------------------------------------------------------------------
  refined <- NULL
  refined_score <- NA_real_

  if (isTRUE(refine)) {
    hs_final <- sort(unique(run_res$hs_seq))
    idx_s <- which.min(abs(hs_final - best$h_s))
    lb_s <- hs_final[max(1, idx_s - 1)]
    ub_s <- hs_final[min(length(hs_final), idx_s + 1)]

    do_p <- use_parallel
    message("Final Refinement (Golden Search)...")

    if (is_spatial_only) {
      refined <- tryCatch(
        golden_search_bandwidth_in(
          formula = formula, Ht = NULL,
          data = HEAVY$data, coords = HEAVY$coords, fixed_vars = HEAVY$fixed_vars,
          kernels = kernels, Model = Model,
          control = control_o,
          ncore = ifelse(do_p, ncore, 1),
          show_progress = show_progress,
          lower.bound = lb_s, upper.bound = ub_s,
          tolerance = tol[1]
        ),
        error = function(e) NULL
      )
    } else {
      ht_final <- sort(unique(run_res$ht_seq))
      idx_t <- which.min(abs(ht_final - best$h_t))
      lb_t <- ht_final[max(1, idx_t - 1)]
      ub_t <- ht_final[min(length(ht_final), idx_t + 1)]

      refined <- tryCatch(
        golden_search_2d_bandwidth(
          formula = formula, data = HEAVY$data, coords = HEAVY$coords, fixed_vars = HEAVY$fixed_vars,
          kernels = kernels, Model = Model, control = control_o,
          lower.bound.space = lb_s, upper.bound.space = ub_s,
          lower.bound.time  = lb_t, upper.bound.time  = ub_t,
          tolerance_s = tol[1], tolerance_t = tol[2]
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(refined) && !is.null(refined$model)) {
      crit <- control_o$criterion %||% "AICc"
      refined_score <- .get_criterion_value(refined$model, crit)
    }
  }

  # --- Compare score vs score (coherent) ---
  refined_ok <- is.finite(refined_score)
  choose_refine <- refined_ok && is.finite(best$score) && (best$score > refined_score)

  best_model <- if (choose_refine) refined$model else {
    H_final <- if (is_spatial_only) best$h_s else c(best$h_s, best$h_t)

    if (!is_spatial_only && best$h_s >= max(run_res$hs_seq) && best$h_t >= max(run_res$ht_seq)) {
      try(MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = NULL, Model = "OLS", control = control_o), silent = TRUE)
    } else {
      try(MGWRSAR(formula, HEAVY$data, HEAVY$coords, HEAVY$fixed_vars, kernels, H = H_final, Model = Model, control = control_o), silent = TRUE)
    }
  }

  elapsed_time <- difftime(Sys.time(), start_time, units = "mins")

  if (verbose) {
    crit_lab <- control$criterion %||% "AICc"
    if (choose_refine) {
      message(sprintf("Chosen refined model: %s=%.6g (improved over grid best %.6g)", crit_lab, refined_score, best$score))
    } else {
      message(sprintf("Chosen grid best model: %s=%.6g", crit_lab, best$score))
    }
  }

  if (acc$isolated_hits > 0) warning(sprintf("\n[WARNING] Isolated points in %d evaluations.\n", acc$isolated_hits))

  list(
    mode =  if (control$Type=='GD') "spatial" else if(control$Type=='T') "temporal" else if(control$Type=='GDT') "spatio-temporal",
    criterion = control$criterion %||% "AICc",
    results = run_res$results,
    best = best,
    refined = refined,
    objective = if (isTRUE(refine)) refined$objective else best$score,
    minimum = if (isTRUE(refine)) refined$minimum else H_final,
    best_model = best_model,
    ctime = elapsed_time
  )
}


