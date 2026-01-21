#' Multiscale Geographically Weighted Regression (MGWR)
#'
#' @description
#' This function estimates a Multiscale Geographically Weighted Regression (MGWR) model
#' based on the proposition of Fotheringham et al. (2017). Unlike standard GWR where a
#' single bandwidth is used for all covariates, MGWR allows for covariate-specific
#' bandwidths. It uses a backfitting algorithm to iteratively estimate the optimal
#' bandwidth and coefficients for each explanatory variable.
#'
#' @usage multiscale_gwr(formula, data, coords, kernels = 'bisq',
#'                       control_mgwr = list(), control = list())
#'
#' @param formula A formula object specifying the model (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing the variables in the model.
#' @param coords A matrix or data frame of coordinates (2 columns for spatial, 3 for spatio-temporal).
#' @param kernels A character string specifying the kernel type.
#' Options include \code{'bisq'} (default), \code{'gauss'}, \code{'triangle'}, \code{'tricube'}, \code{'rectangle'}.
#' @param control_mgwr A named list of control parameters specific to the MGWR backfitting algorithm.
#' See 'Details' for available components.
#' @param control A named list of standard control arguments passed to the internal GWR estimation steps.
#' See 'Details' for available components.
#'
#' @details
#'
#' \strong{Components for \code{control_mgwr}:}
#' \describe{
#'   \item{\code{init}}{Character. The type of model used for initialization. Options are \code{'GWR'} (default) or \code{'lm'} (OLS).}
#'   \item{\code{maxiter}}{Integer. Maximum number of backfitting iterations. Default is 20.}
#'   \item{\code{tolerance}}{Numeric. Convergence threshold based on the change in RMSE or bandwidths. Default is \code{1e-6}.}
#'   \item{\code{nstable}}{Integer. Number of consecutive iterations where bandwidths must remain stable to declare convergence. Default is 6.}
#'   \item{\code{H0}}{Numeric vector. Optional initial bandwidths for each covariate. If \code{NULL}, they are initialized via GWR or global search.}
#'   \item{\code{get_AIC}}{Logical. If \code{TRUE}, calculates the corrected Akaike Information Criterion (AICc) at the end. Default is \code{FALSE}.}
#'   \item{\code{verbose}}{Logical. If \code{TRUE}, prints progress information during backfitting. Default is \code{FALSE}.}
#' }
#'
#' \strong{Components for \code{control}:}
#' \describe{
#'   \item{\code{adaptive}}{Logical. If \code{TRUE} (default), uses an adaptive bandwidth (k-nearest neighbors). If \code{FALSE}, uses a fixed distance bandwidth.}
#'   \item{\code{Type}}{Character. The type of spatial weighting. \code{'GD'} (Geographical Distance, default) or \code{'GDT'} (Geo-Temporal).}
#'   \item{\code{NN}}{Integer. Maximum number of neighbors for matrix truncation (speeds up computation). Default is \code{nrow(data)}.}
#'   \item{\code{ncore}}{Integer. Number of cores to use for parallel computation.}
#'   \item{\code{isgcv}}{Logical. If \code{TRUE}, computes Leave-One-Out Cross-Validation scores. Default is \code{FALSE}.}
#' }
#'
#' @return An object of class \code{mgwrsar} containing:
#' \item{Betav}{Matrix of estimated spatially varying coefficients.}
#' \item{H}{Vector of final optimal bandwidths for each covariate.}
#' \item{RMSE}{Root Mean Square Error of the final model.}
#' \item{residuals}{Vector of residuals.}
#' \item{fitted.values}{Vector of fitted values.}
#' \item{AICc}{Corrected AIC (if \code{get_AIC = TRUE}).}
#' \item{R2}{R-squared of the model.}
#'
#' @references
#' Fotheringham, A. S., Yang, W., & Kang, W. (2017). Multiscale geographically weighted regression (MGWR).
#' \emph{Annals of the American Association of Geographers}, 107(6), 1247-1265.
#'
#' @seealso \code{\link{MGWRSAR}}, \code{\link{TDS_MGWR}}, \code{\link{golden_search_bandwidth}}
#' @export
multiscale_gwr <- function(formula, data, coords, kernels = 'bisq', control_mgwr = list(), control = list()) {

  start <- proc.time()
  n = nrow(data)
  mf <- model.frame(formula, data)
  data <- data[, names(mf)]
  mt <- attr(x = mf, which = "terms")
  Y <- model.extract(mf, "response")
  X = model.matrix(object = mt, data = mf)
  namesX = colnames(X)
  if (colnames(X)[1] == "(Intercept)") colnames(X)[1] <- namesX[1] <- 'Intercept'
  data$Intercept = rep(1, n)
  K = length(namesX)

  # --- control_mgwr (Defaults + Merge) ---
  defaults_mgwr <- list(
    init = 'GWR',
    maxiter = 50,
    nstable = 5,
    tolerance = 0.000001,
    tol=0.001,
    ncore = 1,
    show_progress=FALSE,
    HF = NULL,
    H0 = NULL,
    Ht = NULL,
    Model = NULL,
    model = NULL,
    get_AIC = FALSE,
    verbose = FALSE
  )
  control_mgwr <- modifyList(defaults_mgwr, as.list(control_mgwr))
  list2env(control_mgwr, envir = environment())

  # --- init control ---
  if(is.null(control$NN)) control$NN<-NN<-n else NN<-control$NN
  if(is.null(control$TP)) TP<-control$TP<-1:NN
  controlv <- control
  control$Type = 'GD'
  Model = 'GWR'
  if (get_AIC) {
    control$get_Rk = TRUE
    control$get_ts = TRUE
    control$get_s = TRUE
  }
  controlv$criterion = 'AICc'
  controlv$get_Rk = FALSE
  controlv$get_s = FALSE
  controlv$get_ts = TRUE

  if(verbose) cat("GWR estimation as starting Model \n")


  if (!('indexG' %in% names(control))) {
    while (sum(duplicated(coords)) > 0) {
      set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")
      coords <- jitter(coords, 0.0000001)
    }
    G <- prep_d(coords = coords, NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
    controlv$indexG <- control$indexG <- G$indexG
    controlv$dists <- control$dists <- G$dists
  } else if (!is.null(model)) {
    if (length(model@indexG) == 0) {
      G <- prep_d(coords = coords, NN = control$NN, TP = control$TP, kernels = kernels, Type = control$Type)
      controlv$indexG <- control$indexG <- G$indexG
      controlv$dists <- control$dists <- G$dists
    } else {
      controlv$indexG <- control$indexG <- model@indexG
      controlv$dists <- control$dists <- model@dists
    }
  }

  if (!is.null(model)) {
    HF = model$H
    BETA = model$Betav
    data$eps <- model$residuals
  } else if (init == 'lm') {
    model0 = lm(formula, data)
    BETA = matrix(rep(coef(model0), each = nrow(data)), byrow = FALSE, ncol = length(coef(model0)))
    data$eps <- residuals(model0)
    H = rep(n, K)
  } else {
    if (is.null(H0)) {
      if (control$adaptive) {
        if (kernels[1] == 'gauss') lower.bound = 2 else lower.bound = 2 * K
        upper.bound = control$NN - 2
        tolerance_GWR = 1
      } else {
        imin <- which.min(coords[, 1] + coords[, 2])
        imax <- which.max(coords[, 1] + coords[, 2])
        max_dist = sqrt((coords[imin, 1] - coords[imax, 1])^2 + (coords[imin, 2] - coords[imax, 2])^2)
        lower.bound = 0
        upper.bound = max_dist
        tolerance_GWR = tolerance
      }
      res<- search_bandwidths(formula = formula,
                                data = data,
                                coords = coords,
                                kernels = kernels,
                                Model = Model,
                                control = control,
                                hs_range = c(lower.bound, upper.bound),
                                ht_range =NULL,
                                n_seq = 10,
                                ncore = ncore,
                                n_rounds = 0,
                                tol=0.001,
                                refine=T,
                                show_progress =F,
                                verbose = F
      )
      res$minimum=res$best_model@H
      res$model=res$best_model
      res$objective=res$best_model@AICc
      H0 = res$minimum
      if (get_AIC) {
        modelGWR <- MGWRSAR(formula = formula, H = c(H0, Ht), data = data, coords = coords, fixed_vars = NULL, kernels = kernels, Model = Model, control = control)
      } else  modelGWR <- res$model
    }
    control$get_Rk = FALSE
    BETA = modelGWR@Betav
    H = rep(modelGWR@H, K)
    data$eps <- modelGWR@residuals
    if (verbose) cat('H0=', H0, '\n')
    if (get_AIC) {
      St <- modelGWR@Shat
      Rkk <- Rk <- modelGWR@R_k
      AICc <- aicc_f(modelGWR@residuals, sum(diag(St)), n)
      if (verbose) cat('AICc=', AICc, ' RMSE=', modelGWR@RMSE, '\n')
    }
  }

  drmse <- delta_rmse <- sqrt(mean(Y^2))
  isgcv = FALSE
  iter = 0
  stable <- rep(0, K)
  if (!is.null(HF)) H = HF
  HBETA = list()
  rmse_list <- c()

  # A SUPPRIMER
  # bsl<<-list()
  # A SUPPRIMER

  while ((abs(delta_rmse) > tolerance | any(stable < nstable)) & iter <= maxiter & length(unique(tail(rmse_list))) < 2) {
    iter = iter + 1
    if (verbose) cat('\n')
    if (verbose) cat(' backfitting ')
    for (k in 1:K) {
      var = namesX[k]
      if (verbose)  cat(' ', var, ' ')
      # Utilisation de as.matrix standard au lieu de pipe %>% pour limiter les dépendances si nécessaire
      data$epst <- data$eps + as.matrix(BETA[, k] * data[, var], ncol = 1)
      myformula = as.formula(paste0('epst~', var, '-1'))

      if (is.null(HF)) {
        if (stable[k] < nstable) {
          res<- search_bandwidths(formula = myformula,
                                    data = data,
                                    coords = coords,
                                    kernels = kernels,
                                    Model = Model,
                                    control = controlv,
                                    hs_range = c(lower.bound, upper.bound),
                                    ht_range =NULL,
                                    n_seq = 10,
                                    tol=0.001,
                                    ncore = ncore,
                                    n_rounds =0,
                                    refine=T,
                                    show_progress =F,
                                    verbose = F
          )
          res$minimum=res$best_model@H
          res$model=res$best_model
          res$objective=res$best_model@AICc

          if (get_AIC) {
            res$model <- MGWRSAR(formula = myformula, H = c(res$minimum, Ht), data = data, coords = coords, fixed_vars = NULL, kernels = kernels, Model = Model, control = control)
          }
          modellm <- MGWRSAR(formula = myformula, data = data, coords = coords, fixed_vars = NULL, kernels = kernels, H = H0, Model = 'OLS', control = control)
          AIClm = AIC(lm(myformula, data))
          if (res$model@RMSE > modellm@RMSE) {
            res$minimum = n
            res$objective = AIClm
          }
          H[k] <- h <- res$minimum
        } else h = H[k]
      } else h = HF[k]

      if (h == n) {
        if (stable[k] >= nstable) {
          modellm <- MGWRSAR(formula = myformula, data = data, coords = coords, fixed_vars = NULL, kernels = kernels, H = H0, Model = 'OLS', control = control)
          AIClm = AIC(lm(myformula, data))
        }
        modelGWR = modellm
        modelGWR@H = c(n, Ht)
        modelGWR@Betav = as.matrix(rep(modelGWR@Betac, n), ncol = 1)
        modelGWR@residuals = residuals(modellm)
      } else {
        if (stable[k] < nstable & is.null(HF)) modelGWR <- res$model else modelGWR <- MGWRSAR(formula = myformula, data = data, coords = coords, fixed_vars = NULL, kernels = kernels, H = c(h, Ht), Model = 'multiscale_gwr', control = control)
      }
      BETA[, k] <- modelGWR@Betav
      data$eps = modelGWR@residuals
      if (get_AIC) {
        Sk <- modelGWR@Shat
        Rkk[[k]] <- eigenMapMatMult(Sk, Rk[[k]]) + Sk - eigenMapMatMult(Sk, St)
        St = St + Rkk[[k]] - Rk[[k]]
        new_ts = sum(diag(St))
        AICg <- aicc_f(data$eps, new_ts, n)
        Rk[[k]] <- Rkk[[k]]
      }
      if (iter > 1) if (oldH[k] == H[k]) stable[k] = stable[k] + 1 else stable[k] = 0
    }
    HBETA[[iter]] <- BETA
    rmse_list <- c(rmse_list, rmse)
    oldH <- H
    delta_rmse = (drmse - sqrt(mean(data$eps^2))) / drmse
    drmse = sqrt(mean(data$eps^2))
    if (verbose & get_AIC) cat('\n iter ', iter, ' H ', H, ' RMSE=', sqrt(mean(data$eps^2)), ' delta_rmse ', delta_rmse, ' AICc=', AICg, ' stable', stable, '\n') else  if (verbose) cat('\n iter ', iter, ' H ', H, ' RMSE=', sqrt(mean(data$eps^2)), ' delta_rmse ', delta_rmse, ' stable', stable, '\n')
  }
  if (verbose) cat('Time = ', (proc.time() - start)[3], '\n')
  modelGWR@Model = 'multiscale_gwr'
  modelGWR@Betav = BETA
  if (get_AIC) modelGWR@AICc = AICg
  modelGWR@residuals = data$eps
  modelGWR@fit = Y - data$eps
  modelGWR@RMSE = sqrt(mean(modelGWR@residuals^2))
  modelGWR@RMSEtp = sqrt(mean(modelGWR@residuals^2))
  modelGWR@H = c(H, Ht)
  modelGWR@X = X
  modelGWR@ctime = (proc.time() - start)[3]
  modelGWR
}
