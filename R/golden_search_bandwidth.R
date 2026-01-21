#' @title Golden search bandwidth (deprecated)
#' @description
#' `golden_search_bandwidth()` is deprecated. Please use [search_bandwidths()],
#' which provides a unified interface for spatial and spatio-temporal bandwidth
#' search and can run golden-search refinement internally.
#'
#' This wrapper preserves backward compatibility by translating the historical
#' arguments (`lower.bound`, `upper.bound`, `Ht`) to the new interface.
#'
#' @inheritParams search_bandwidths
#' @param Ht Optional temporal bandwidth (used only for legacy calls). If not
#'   `NULL`, it is mapped to `ht_range = c(Ht, Ht)` (i.e., fixed temporal
#'   bandwidth) when `n_rounds = 0`.
#' @param lower.bound,upper.bound Numeric bounds for the spatial bandwidth
#'   search; mapped to `hs_range`.
#' @param tolerance Numeric tolerance; mapped to `tol` when relevant.
#' @param ncore Number of cores to use.
#' @param show_progress Logical; forwarded to [search_bandwidths()].
#'
#' @return A list returned by [search_bandwidths()] (see its Value section).
#' @export
golden_search_bandwidth <- function(
    formula, Ht = NULL, data, coords, fixed_vars = NULL, kernels,
    Model = "GWR", control, lower.bound, upper.bound,
    tolerance = 1e-6,
    ncore = 1,
    show_progress = TRUE
) {
  # --- Deprecation message (no hard dependency) ---
  if (requireNamespace("lifecycle", quietly = TRUE)) {
    lifecycle::deprecate_warn(
      "1.3.#", # <- mets ta version de dépréciation
      "golden_search_bandwidth()",
      "search_bandwidths()"
    )
  } else {
    warning("`golden_search_bandwidth()` is deprecated; use `search_bandwidths()` instead.",
            call. = FALSE)
  }

  ptm <- proc.time()
  `%||%` <- function(x, y) if (!is.null(x)) x else y

  get_SFdata()

  mf <- model.frame(formula, data)
  data <- data[,names(mf)]

  acc <- new.env(parent = emptyenv())
  acc$isolated_hits <- 0L

  # -- control setup
  control_o <- control
  control$SE      <- FALSE
  control$get_s   <- FALSE
  control$get_Rk  <- FALSE
  if(is.null(control$Type)) control$Type <- 'GD'

  # -- jitter coords
  if(control$Type %in% c('GD','GDT')) coords<-make_unique_by_structure(coords)
  if(control$Type %in% c('GDT','T')) control$Z<-make_unique_by_structure(control$Z)

  adaptive <- isTRUE(control$adaptive[1])
  golden.ratio <- 2 / (sqrt(5) + 1)

  # -- starting points
  x1 <- upper.bound - golden.ratio * (upper.bound - lower.bound)
  x2 <- lower.bound + golden.ratio * (upper.bound - lower.bound)
  if (adaptive) { x1 <- floor(x1); x2 <- ceiling(x2) }

  # -- distances computation (Lourd !)
  if (is.null(control$dists)) {
    if (is.null(control$NN)) control$NN <- nrow(data)
    if (is.null(control$TP)) control$TP <- 1:nrow(data)
    S <- if (control$Type == 'T') as.matrix(control$Z, ncol=1)
    else if(control$Type == 'GDT') as.matrix(cbind(coords, control$Z), ncol=3)
    else as.matrix(coords)

    stage1 <- prep_d(coords = S, NN = control$NN, TP = control$TP,
                     kernels = kernels, Type = control$Type)
    control$indexG <- stage1$indexG
    control$dists  <- stage1$dists
  }
  control_o <- control

  # =====================================================
  # safe_eval definition
  # =====================================================
  safe_eval <- function(H) {
    crit <- control$criterion %||% "AICc"

    # ---  CV isgcv = TRUE ---
    if (crit %in% c("CV", "CVtp")) {
      control2 <- control
      control2$isgcv <- TRUE

      mod <- tryCatch(
        MGWRSAR(formula = formula, data = data, coords = coords, fixed_vars = fixed_vars,
                kernels = kernels, H = H, Model = Model, control = control2),
        error = function(e) NULL
      )
      if (is.null(mod)) return(Inf)
      slots <- slotNames(mod)

      if (crit == "CV" && "RMSE" %in% slots) val <- mod@RMSE
      else if (crit == "CVtp" && "RMSEtp" %in% slots) val <- mod@RMSEtp
      else val <- Inf

    } else {
      control_AIC <- control
      control_AIC$get_ts <- TRUE

      mod <- tryCatch(
        MGWRSAR(formula = formula, data = data, coords = coords, fixed_vars = fixed_vars,
                kernels = kernels, H = H, Model = Model, control = control_AIC),
        error = function(e) NULL
      )

      if (is.null(mod)) return(Inf)
      slots <- slotNames(mod)

      if (crit %in% slots) val <- slot(mod, crit)
      else if ("AICc" %in% slots) val <- mod@AICc
      else if ("AIC" %in% slots) val <- mod@AIC
      else val <- Inf
    }
    return(as.numeric(val))
  }

  # =====================================================
  # OPTIMISATION LINUX vs MAC/WIN
  # =====================================================
  cl <- NULL
  use_parallel <- (ncore > 1) && requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)

  if (use_parallel) {
    sysname <- Sys.info()[['sysname']]

    # Cas 1 : Linux (Forking) - Très rapide, mémoire partagée
    if (sysname == "Linux") {
      doParallel::registerDoParallel(cores = ncore)
      # Pas besoin d'exporter les données, le fork hérite de l'environnement !

    } else {
      # Cas 2 : Windows/Mac (Sockets) - Plus lent, nécessite copie
      cl <- parallel::makeCluster(ncore)
      doParallel::registerDoParallel(cl)

      # Export explicite uniquement pour les Sockets
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages(library(mgwrsar))
        NULL
      })
      # On exporte tout ce dont safe_eval a besoin
      # Note: 'control' est lourd mais exporté une seule fois ici
      parallel::clusterExport(cl, varlist = c("formula", "data", "coords", "fixed_vars",
                                              "kernels", "Model", "control", "safe_eval"),
                              envir = environment())
    }

    # Nettoyage automatique en fin de fonction
    on.exit({
      if (!is.null(cl)) parallel::stopCluster(cl)
      doParallel::stopImplicitCluster()
    }, add = TRUE)
  }

  # Fonction d'évaluation qui utilise le backend enregistré (Fork ou Socket)
  eval_pair_parallel <- function(h1, h2) {
    H_list <- list(c(h1, Ht), c(h2, Ht))

    if (use_parallel) {
      # Utilisation de foreach qui profite du backend enregistré
      # .export est NULL car géré soit par Fork (Linux), soit par clusterExport (Mac/Win)
      res_list <- foreach::foreach(h = H_list, .combine = c,
                                   .packages = c('mgwrsar')) %dopar% {
                                     safe_eval(h)
                                   }
    } else {
      res_list <- sapply(H_list, safe_eval)
    }
    return(res_list)
  }

  # -- first evaluations
  f12 <- eval_pair_parallel(x1, x2)
  f1 <- f12[1]; f2 <- f12[2]

  # ========================
  # golden search loop
  # ========================
  iteration <- 0L
  if (adaptive) tolerance <- 1
  x1_p <- 1; x2_p <- 2
  max_iter <- 50L

  if(show_progress) {
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  while (abs(upper.bound - lower.bound) > tolerance &&
         ((adaptive && abs(x1 - x2) > 1 && (x1_p != x1 || x2_p != x2)) || !adaptive)) {

    iteration <- iteration + 1L
    if(show_progress) setTxtProgressBar(pb, iteration)
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

  # -- Final Result
  res <- (lower.bound + upper.bound) / 2
  if (adaptive) {
    candidates <- unique(c(x1, x2, lower.bound, upper.bound))
    # Pour 4 candidats max, pas besoin de paralélisme (trop d'overhead)
    scores <- vapply(candidates, function(h) safe_eval(c(h, Ht)), numeric(1))
    res <- candidates[which.min(scores)]
    objective <- min(scores)
  } else {
    objective <- safe_eval(c(res, Ht))
  }

  # -- final model
  final_model <- MGWRSAR(
    formula, data, coords, fixed_vars, kernels,
    H = c(res, Ht), Model = Model, control = control_o
  )

  # -- warning handling
  if (acc$isolated_hits > 0) {
    warning(sprintf(
      "\n[WARNING] Non-adaptive %s kernel produced isolated points in %d evaluations.\n[NOTE] Consider using adaptive or Gaussian kernels.",
      kernels[1], acc$isolated_hits
    ))
  }

  elapsed_time <- (proc.time() - ptm)[3]

  list(
    mode = if (control$Type=='GD') "spatial" else if(control$Type=='T') "temporal" else if(control$Type=='GDT') "spatio-temporal",
    minimum=res,
    criterion = control$criterion %||% "AICc",
    objective = objective,
    best_model = final_model,
    ctime = elapsed_time
  )
}


