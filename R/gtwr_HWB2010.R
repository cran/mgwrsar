#' Control keys consumed by gtwr_HWB2010 itself
#'
#' These keys are removed from \code{control} before it is forwarded to
#' \code{\link{MGWRSAR}}, which warns about unknown control names.
#' @noRd
.gtwr_own_keys <- c("causal", "tail_tol", "t.units")

#' Keys of control that gtwr_HWB2010 sets itself and refuses to have overridden
#' @noRd
.gtwr_locked_keys <- c("Type", "kernels", "adaptive", "alpha", "Z")


#' Convert Huang's (h_ST, tau) into the mgwrsar (h_S, h_T) pair
#'
#' Huang et al. (2010) weight observation \eqn{j} around \eqn{i} with
#' \eqn{\alpha_{ij} = \exp\{-[d^2_S + \tau d^2_T]/h^2_{ST}\}} (their equations
#' 11-12, with the scale convention \eqn{\lambda = 1}, \eqn{\tau = \mu/\lambda}).
#' The gaussian kernel of mgwrsar is \eqn{\exp\{-\frac{1}{2}(d/h)^2\}}, hence the
#' factor \eqn{\sqrt{2}}. Matching the two expressions term by term gives
#' \eqn{h_S = h_{ST}/\sqrt{2}} and \eqn{h_T = h_S/\sqrt{\tau}}.
#'
#' The two degenerate cases are the special cases of Huang's equation 10 rather
#' than limits of the \eqn{\tau} parameterisation: \code{tau = 0} is
#' \eqn{\mu = 0} (pure GWR, flat temporal kernel) and \code{tau = Inf} is
#' \eqn{\lambda = 0} (pure TWR, flat spatial kernel), in which case
#' \code{h_st} is read on the temporal axis.
#'
#' Use it to turn a \eqn{(h_{ST}, \tau)} pair into the \code{hs_range} and
#' \code{ht_range} that \code{\link{search_bandwidths}} expects, and
#' \code{\link{bw_gdt2hwb}} to read its result back.
#'
#' @param h_st numeric scalar, Huang's single spatio-temporal bandwidth.
#' @param tau numeric scalar in \code{[0, Inf]}, Huang's ratio \eqn{\mu/\lambda}.
#' @return A numeric vector \code{c(h_S, h_T)}.
#' @seealso bw_gdt2hwb, gtwr_HWB2010, search_bandwidths
#' @export
#' @examples
#' bw_hwb2gdt(h_st = 0.0849, tau = 0.0225)   # -> c(0.06, 0.4)
bw_hwb2gdt <- function(h_st, tau) {
  h <- h_st / sqrt(2)
  if (is.infinite(tau)) return(c(Inf, h))
  if (tau == 0)         return(c(h, Inf))
  c(h, h / sqrt(tau))
}

#' Convert an mgwrsar (h_S, h_T) pair back into Huang's (h_ST, tau)
#'
#' Inverse of \code{\link{bw_hwb2gdt}}: reads the optimum returned by
#' \code{\link{search_bandwidths}} in the parameterisation of
#' Huang et al. (2010), with \eqn{h_{ST} = \sqrt{2} h_S} and
#' \eqn{\tau = (h_S/h_T)^2}.
#'
#' @param h_s numeric scalar, spatial bandwidth on the mgwrsar scale.
#' @param h_t numeric scalar, temporal bandwidth on the mgwrsar scale.
#' @return A numeric vector \code{c(h_st, tau)}.
#' @seealso bw_hwb2gdt, as_gtwr, search_bandwidths
#' @export
#' @examples
#' bw_gdt2hwb(h_s = 0.06, h_t = 0.4)         # -> c(0.0849, 0.0225)
bw_gdt2hwb <- function(h_s, h_t) {
  if (is.infinite(h_s)) return(c(sqrt(2) * h_t, Inf))
  c(sqrt(2) * h_s, if (is.infinite(h_t)) 0 else (h_s / h_t)^2)
}


#' Largest kernel weight carried by the outer ring of the kNN screening
#'
#' \code{prep_d} pre-selects \code{NN} neighbours per target point using a joint
#' space-time metric standardised by sample standard deviations, i.e. with an
#' implicit ratio that is generally not the calibrated \code{tau}. Since the
#' gaussian kernel has unbounded support, a mismatch can truncate the tail of the
#' weight distribution. The focal point carries a weight of exactly 1, so the raw
#' weight of a boundary neighbour is also its weight relative to the row maximum.
#'
#' @param dists the \code{dists} list returned by \code{prep_d}.
#' @param h_s,h_t bandwidths on the mgwrsar scale.
#' @param ring fraction of the retained neighbours treated as the outer ring.
#' @return A numeric scalar in \code{[0, 1]}.
#' @noRd
.gtwr_tail_weight <- function(dists, h_s, h_t, ring = 0.02) {
  dS <- dists[["dist_s"]]
  dT <- dists[["dist_t"]]
  if (is.null(dS) || is.null(dT)) return(NA_real_)
  k <- ncol(dS)
  m <- max(1L, ceiling(ring * k))
  cols <- seq.int(k - m + 1L, k)
  w <- exp(-0.5 * (dS[, cols, drop = FALSE] / h_s)^2 -
           0.5 * (dT[, cols, drop = FALSE] / h_t)^2)
  suppressWarnings(max(w, na.rm = TRUE))
}


#' Geographically and Temporally Weighted Regression (Huang, Wu and Barry 2010)
#'
#' Calibrates the GTWR model of Huang, Wu and Barry (2010) in its own
#' parameterisation: a single spatio-temporal bandwidth \code{h_st} and a
#' scale ratio \code{tau}. The estimator is the gaussian, fixed-bandwidth
#' special case of the space-time kernel of \code{\link{MGWRSAR}}
#' (\code{Type = 'GDT'}), reached through the exact change of variables
#' \eqn{h_S = h_{ST}/\sqrt{2}}, \eqn{h_T = h_S/\sqrt{\tau}}.
#'
#' @param formula a formula.
#' @param data a data.frame.
#' @param coords a n x 2 matrix or data.frame of spatial coordinates. The time
#'   index is passed through \code{time}, never as a third column.
#' @param time a numeric vector of length n holding the time index. Objects of
#'   class \code{Date}, \code{POSIXct} or \code{POSIXlt} are accepted only when
#'   \code{control$t.units} is given (any unit understood by
#'   \code{\link[base]{difftime}}); they are then converted to elapsed time from
#'   the earliest observation, expressed in that unit.
#' @param h_st the spatio-temporal bandwidth of Huang et al. (2010), on the scale
#'   of the combined distance \eqn{d_{ST}}.
#' @param tau the ratio \eqn{\mu/\lambda} of Huang et al. (2010), with the
#'   convention \eqn{\lambda = 1}. \code{tau = 0} gives a plain GWR and
#'   \code{tau = Inf} a plain TWR.
#' @param fixed_vars a character vector of covariates with a constant
#'   coefficient. \code{NULL} (default) gives the model of Huang et al. (2010);
#'   anything else gives a mixed GTWR, which is an extension of this package and
#'   not part of HWB2010.
#' @param control list of extra control arguments forwarded to
#'   \code{\link{MGWRSAR}}, plus the keys described below.
#'
#' @details
#' Three control keys are consumed by this function and are not forwarded:
#' \describe{
#'   \item{causal}{Logical, default \code{FALSE}. \code{FALSE} uses a symmetric
#'     temporal kernel, which is the model of Huang et al. (2010). \code{TRUE}
#'     switches to a past-only temporal kernel (\code{'gauss_past'}), an
#'     extension in the spirit of Wu et al. (2014), not HWB2010.}
#'   \item{tail_tol}{Numeric, default \code{0.01}. Threshold above which a
#'     warning reports that the kNN pre-selection may be truncating the tail of
#'     the gaussian kernel. Only checked when \code{control$NN < n}.}
#'   \item{t.units}{Character, required when \code{time} is a date or time
#'     object, ignored otherwise.}
#' }
#'
#' Setting \code{control$NN} below \code{n} activates the kNN pre-selection of
#' \code{prep_d}, which ranks neighbours in a joint space-time metric
#' standardised by sample standard deviations. That ranking does not depend on
#' \code{tau}, so with screening on the retained neighbourhood is not the one
#' implied by the calibrated \code{tau}: in particular \code{tau = 0} then no
#' longer reproduces a plain GWR exactly. Since the gaussian kernel has no
#' compact support, \code{NN = n} (the default) is the faithful setting, and
#' \code{w_tail} is reported so that any screening can be checked.
#'
#' The keys \code{Type}, \code{kernels}, \code{adaptive}, \code{alpha} and
#' \code{Z} are set by this function and supplying them raises an error. In
#' particular the kernel is locked to gaussian, because the separability
#' \eqn{W^{ST} = W^S \odot W^T} that Huang's equation 12 relies on holds for no
#' other kernel, and adaptive bandwidths are refused because mgwrsar adapts each
#' dimension separately whereas Huang ranks neighbours in the joint space-time
#' metric. Use \code{MGWRSAR(Type = 'GDT')} directly for those variants.
#'
#' Bandwidths are not calibrated here. Use \code{\link{search_bandwidths}},
#' the calibration entry point of the package for \code{Type = 'GDT'}, with
#' \code{hs_range}/\code{ht_range} obtained from \code{\link{bw_hwb2gdt}},
#' then read its optimum back with \code{\link{bw_gdt2hwb}} or turn its
#' \code{best_model} into a \code{gtwr} with \code{\link{as_gtwr}}.
#'
#' @return An object of class \code{\link{gtwr-class}}, which extends
#'   \code{\link{mgwrsar-class}} and therefore supports \code{summary},
#'   \code{coef}, \code{fitted}, \code{residuals}, \code{plot} and
#'   \code{predict}. Out-of-sample prediction requires
#'   \code{method_pred} in \code{c('model','tWtp_model')} and
#'   \code{newdata_coords} with columns (x, y, t).
#'
#' @references
#' Huang, B., Wu, B. and Barry, M. (2010). Geographically and temporally weighted
#' regression for modeling spatio-temporal variation in house prices.
#' International Journal of Geographical Information Science, 24(3), 383-401.
#'
#' Wu, B., Li, R. and Huang, B. (2014). A geographically and temporally weighted
#' autoregressive model with application to housing prices. International Journal
#' of Geographical Information Science, 28(6), 1186-1204.
#'
#' @seealso MGWRSAR, search_bandwidths, as_gtwr, bw_hwb2gdt, predict.mgwrsar
#' @export
#' @examples
#' \donttest{
#' library(mgwrsar)
#' data(mydata)
#' coords <- as.matrix(mydata[, c("x", "y")])
#' mydata$time <- rep(1:10, length.out = nrow(mydata))
#'
#' ## coords are in metres here, hence a bandwidth of a few kilometres;
#' ## tau rescales the time axis so that one period weighs like h_st/sqrt(tau)
#' m <- gtwr_HWB2010(formula = 'Y_gwr ~ X1 + X2 + X3', data = mydata,
#'                   coords = coords, time = mydata$time,
#'                   h_st = 4000, tau = 1e6,
#'                   control = list(SE = TRUE))
#' summary(m)
#'
#' ## out-of-sample prediction on (x, y, t)
#' newc <- cbind(coords[1:10, ], mydata$time[1:10])
#' predict(m, newdata = mydata[1:10, ], newdata_coords = newc,
#'         method_pred = 'model')
#' }
gtwr_HWB2010 <- function(formula, data, coords, time, h_st, tau,
                         fixed_vars = NULL, control = list()) {

  mycall <- match.call()

  if (missing(h_st) || length(h_st) != 1 || !is.finite(h_st) || h_st <= 0)
    stop("h_st must be a single finite positive value")
  if (missing(tau) || length(tau) != 1 || is.na(tau) || tau < 0)
    stop("tau must be a single value in [0, Inf]")

  st <- .gtwr_prepare(data, coords, time, control, caller = "gtwr_HWB2010")

  H <- bw_hwb2gdt(h_st, tau)
  w_tail <- .gtwr_tail_weight(st$control$dists, H[1], H[2])
  if (st$NN < st$n && is.finite(w_tail) && w_tail > st$tail_tol)
    warning("kNN pre-selection may truncate the gaussian tail: the outer ring ",
            "of the NN = ", st$NN, " retained neighbours still carries a weight of ",
            signif(w_tail, 3), " (tail_tol = ", st$tail_tol,
            "). Increase control$NN.", call. = FALSE)

  res <- MGWRSAR(formula = formula, data = data, coords = st$coords,
                 fixed_vars = fixed_vars, kernels = st$kernels, H = H,
                 Model = if (is.null(fixed_vars)) "GWR" else "MGWR",
                 control = st$control)

  ## Slots are copied with check = FALSE: an mgwrsar object left by MGWRSAR
  ## can hold the empty placeholder of the virtual "Matrix" slot W, which does
  ## not pass validation when re-assigned.
  out <- new("gtwr")
  for (s in methods::slotNames("mgwrsar"))
    methods::slot(out, s, check = FALSE) <- methods::slot(res, s)

  out@h_st   <- h_st
  out@tau    <- tau
  out@causal <- st$causal
  out@w_tail <- as.numeric(w_tail)
  out@mycall <- mycall
  invisible(out)
}


#' Shared validation and control assembly for the GTWR wrappers
#'
#' Performs the argument checks, the de-duplication that \code{MGWRSAR} would
#' perform anyway, the locking of the control keys, and the single
#' \code{prep_d} call whose result is reused by every model fit (the kNN
#' pre-selection does not depend on the bandwidths).
#'
#' @param data,coords,time as in \code{\link{gtwr_HWB2010}}.
#' @param control the user control list.
#' @param caller name used in the error messages.
#' @return A list with the prepared \code{coords}, \code{time}, \code{kernels},
#'   \code{control}, and the scalars \code{causal}, \code{tail_tol}, \code{n},
#'   \code{NN}.
#' @noRd
.gtwr_prepare <- function(data, coords, time, control, caller) {

  if (!is.list(control)) stop("control must be a list")
  locked <- intersect(names(control), .gtwr_locked_keys)
  if (length(locked))
    stop(caller, " sets these control keys itself, remove them: ",
         paste(locked, collapse = ", "),
         ". Use MGWRSAR(Type='GDT') if you need to set them.")

  causal   <- isTRUE(control$causal)
  tail_tol <- if (is.null(control$tail_tol)) 0.01 else control$tail_tol
  t.units  <- control$t.units
  control  <- control[setdiff(names(control), .gtwr_own_keys)]

  coords <- as.matrix(coords)
  if (ncol(coords) != 2)
    stop("coords must have exactly 2 columns; the time index goes in `time`, ",
         "not in a third column of coords")
  if (!is.numeric(coords)) stop("coords must be numeric")

  n <- nrow(data)
  if (nrow(coords) != n)
    stop("coords and data must have the same number of rows")

  if (missing(time)) stop("time must be provided")
  if (inherits(time, c("Date", "POSIXct", "POSIXlt"))) {
    if (is.null(t.units))
      stop("time is a date/time object: provide control$t.units ",
           "(e.g. 'days', 'weeks') to fix the unit of the temporal distance")
    time <- as.numeric(difftime(time, min(time), units = t.units))
  }
  time <- as.numeric(time)
  if (length(time) != n)
    stop("time must have the same length as the number of rows of data")
  if (anyNA(time)) stop("time must not contain NA")
  if (anyNA(coords)) stop("coords must not contain NA")

  kernels <- c("gauss", if (causal) "gauss_past" else "gauss")

  ## ------------------------------------------------------------------
  ## same de-duplication as MGWRSAR, applied once so the pre-computed
  ## distances and the ones MGWRSAR would build cannot diverge
  ## (make_unique_by_structure is idempotent)
  ## ------------------------------------------------------------------
  coords <- make_unique_by_structure(coords)
  time   <- make_unique_by_structure(time)

  control$Type     <- "GDT"
  control$Z        <- time
  control$adaptive <- c(FALSE, FALSE)
  control$alpha    <- 1

  NN <- if (is.null(control$NN)) n else min(as.integer(control$NN), n)
  control$NN <- NN
  ## TP is pinned down here so that the pre-computed distances and the control
  ## handed to MGWRSAR or to golden_search_2d_bandwidth describe the same points
  TP <- if (is.null(control$TP)) seq_len(n) else control$TP
  control$TP <- TP
  extrapol <- isTRUE(control$TP_estim_as_extrapol)

  ## ------------------------------------------------------------------
  ## distances computed once, reused by MGWRSAR and by the tail check
  ## ------------------------------------------------------------------
  if (is.null(control$dists) || is.null(control$indexG)) {
    stage1 <- prep_d(coords = as.matrix(cbind(coords, time)),
                     NN = NN, TP = TP, extrapol = extrapol,
                     kernels = kernels, Type = "GDT")
    control$dists  <- stage1$dists
    control$indexG <- stage1$indexG
  }

  list(coords = coords, time = time, kernels = kernels, control = control,
       causal = causal, tail_tol = tail_tol, n = n, NN = NN)
}


#' Coerce a calibrated Type='GDT' model into Huang's parameterisation
#'
#' Bandwidth calibration for GTWR goes through \code{\link{search_bandwidths}},
#' the entry point of the package for space-time kernels, which optimises on the
#' \code{(h_S, h_T)} axes. This turns its \code{best_model} into an object of
#' class \code{\link{gtwr-class}}, so that the result reads in the
#' \eqn{(h_{ST}, \tau)} parameterisation of Huang et al. (2010) and inherits
#' the dedicated \code{summary}.
#'
#' The model is checked against the HWB2010 contract: \code{Type = 'GDT'}, a
#' gaussian kernel on both axes, fixed (non-adaptive) bandwidths and a purely
#' multiplicative kernel (\code{alpha = 1}). Anything else is rejected, because
#' the separability that Huang's equation 12 relies on no longer holds.
#'
#' @param object a model of class \code{\link{mgwrsar-class}} fitted with
#'   \code{Type = 'GDT'}, typically \code{search_bandwidths(...)$best_model}.
#' @return An object of class \code{\link{gtwr-class}}. The \code{w_tail} slot
#'   is \code{NA} for a coerced model: the kNN screening diagnostic needs the
#'   distance matrices, which the fitted object does not carry. Refit with
#'   \code{\link{gtwr_HWB2010}} if you need it.
#'
#' @seealso gtwr_HWB2010, search_bandwidths, bw_gdt2hwb
#' @export
#' @examples
#' \donttest{
#' library(mgwrsar)
#' simu <- simu_multiscale(n = 1000, myseed = 1, type = 'GG2024',
#'                         config_beta = 'spatiotemp_old', config_snr = 0.9)
#' dat <- simu$mydata
#'
#' res <- search_bandwidths(
#'   formula = as.formula('Y ~ X1 + X2 + X3'), data = dat,
#'   coords = simu$coords, kernels = c('gauss', 'gauss'), Model = 'GWR',
#'   control = list(Z = dat$time, criterion = 'AICc', Type = 'GDT',
#'                  adaptive = c(FALSE, FALSE), alpha = 1),
#'   hs_range = c(0, 1.4), ht_range = c(0, 1),
#'   n_seq = 12, n_rounds = 3, refine = FALSE)
#'
#' m <- as_gtwr(res$best_model)
#' summary(m)
#' bw_gdt2hwb(res$minimum[1], res$minimum[2])
#' }
as_gtwr <- function(object) {

  if (!is(object, "mgwrsar"))
    stop("object must be of class mgwrsar")
  if (is(object, "gtwr")) return(object)

  if (!identical(object@Type, "GDT"))
    stop("as_gtwr expects a model fitted with Type = 'GDT', got Type = '",
         object@Type, "'")

  k <- object@kernels
  if (length(k) != 2 || k[1] != "gauss" ||
      unlist(strsplit(k[2], "_"))[1] != "gauss")
    stop("HWB2010 requires a gaussian kernel on both axes, got: ",
         paste(k, collapse = ", "),
         ". The separability W^ST = W^S * W^T holds for no other kernel.")
  if (any(object@adaptive))
    stop("HWB2010 requires fixed bandwidths; this model uses an adaptive one ",
         "on at least one axis, which ranks neighbours per dimension rather ",
         "than in the joint space-time metric.")
  if (length(object@alpha) && !isTRUE(all.equal(object@alpha, 1)))
    stop("HWB2010 requires a purely multiplicative kernel (alpha = 1), got ",
         object@alpha)

  hw <- bw_gdt2hwb(object@H[1], object@Ht[1])

  out <- new("gtwr")
  for (s in methods::slotNames("mgwrsar"))
    methods::slot(out, s, check = FALSE) <- methods::slot(object, s)
  out@h_st   <- hw[1]
  out@tau    <- hw[2]
  out@causal <- length(unlist(strsplit(k[2], "_"))) > 1 &&
                unlist(strsplit(k[2], "_"))[2] == "past"
  out@w_tail <- NA_real_
  out
}
