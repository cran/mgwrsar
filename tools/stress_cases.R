# Stress campaign harness: estimator configurations x adversarial generators
# (tools/stress_generators.R). stress_run() fits every pair under tryCatch and
# returns one row per pair: outcome (ok / error / warning text), timing and the
# invariant checks of stress_invariants(). Used interactively to triage, then
# the configurations kept become tests (tests/testthat/test-stress_*.R).

if (!exists("stress_generators")) {
  f <- if (file.exists("tools/stress_generators.R")) "tools/stress_generators.R" else "../../tools/stress_generators.R"
  source(f)
}

# Invariants a fitted mgwrsar object must satisfy whatever the data
stress_invariants <- function(m, d) {
  n <- nrow(d); out <- character()
  chk <- function(ok, msg) if (!isTRUE(ok)) out <<- c(out, msg)
  if (methods::is(m, "mgwrsar")) {
    B <- as.matrix(m@Betav)
    if (length(B)) chk(all(is.finite(B)), "Betav not finite")
    if (length(m@Betac)) chk(all(is.finite(m@Betac)), "Betac not finite")
    if (length(m@fit)) {
      chk(length(m@fit) == n, "fit length != n")
      chk(all(is.finite(m@fit)), "fit not finite")
    }
    if (length(m@residuals) && length(m@fit))
      chk(isTRUE(all.equal(as.numeric(m@residuals), as.numeric(d$Y - m@fit), tolerance = 1e-8)), "residuals != Y - fit")
    if (length(m@TS)) {
      chk(all(is.finite(m@TS)), "TS not finite")
      chk(all(m@TS >= -1e-8 & m@TS <= 1 + 1e-8), sprintf("TS out of [0,1] (min %.3g, max %.3g)", min(m@TS), max(m@TS)))
    }
    if (length(m@tS)) chk(is.finite(m@tS) && m@tS >= 0 && m@tS <= n + 1e-6, sprintf("tS out of [0,n] (%.3g)", m@tS))
    if (length(m@AICc)) chk(is.finite(m@AICc) || (length(m@tS) && n - 1 - m@tS <= 0), "AICc not finite although not saturated")
    if (length(m@sev)) chk(all(is.finite(m@sev) & m@sev >= 0), "sev not finite or negative")
    if (length(m@se)) chk(all(is.finite(m@se) & m@se >= 0), "se not finite or negative")
    # Inf = global bandwidth; a global OLS model (e.g. the corner of a bandwidth
    # search) has no bandwidth at all
    if (!identical(m@Model, "OLS")) {
      if (length(m@H)) chk(all(!is.nan(m@H) & m@H > 0), "H NaN or <= 0")
      if (length(m@Ht)) chk(all(!is.nan(m@Ht) & m@Ht > 0), "Ht NaN or <= 0")
    }
    if (length(m@Shat)) {
      S <- as.matrix(m@Shat)
      if (nrow(S) == n && ncol(S) == n && length(m@TS) == n)
        chk(isTRUE(all.equal(unname(diag(S)), as.numeric(m@TS), tolerance = 1e-8)), "diag(Shat) != TS")
    }
  }
  out
}

# One configuration: a label, the estimator call (a function of the data), and
# a predicate on the data saying whether the configuration applies
stress_config <- function(label, fit, applies = function(d) TRUE, family = "fit")
  list(label = label, fit = fit, applies = applies, family = family)

stress_run <- function(configs, generators = stress_generators(), verbose = TRUE,
                       only_configs = NULL, only_generators = NULL) {
  rows <- list()
  for (gn in names(generators)) {
    if (!is.null(only_generators) && !(gn %in% only_generators)) next
    d <- generators[[gn]]()
    for (cf in configs) {
      if (!is.null(only_configs) && !(cf$label %in% only_configs)) next
      if (!isTRUE(cf$applies(d))) next
      t0 <- Sys.time(); warns <- character()
      res <- withCallingHandlers(
        tryCatch(cf$fit(d), error = function(e) e),
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
      secs <- as.numeric(Sys.time() - t0, units = "secs")
      if (inherits(res, "error")) {
        row <- data.frame(generator = gn, config = cf$label, outcome = "error",
                          detail = substr(conditionMessage(res), 1, 160), secs = secs, stringsAsFactors = FALSE)
      } else {
        inv <- if (methods::is(res, "mgwrsar")) stress_invariants(res, d) else
          if (is.list(res) && !is.null(res$invariants)) res$invariants else character()
        row <- data.frame(generator = gn, config = cf$label,
                          outcome = if (length(inv)) "invariant" else "ok",
                          detail = substr(paste(unique(c(inv, warns)), collapse = " | "), 1, 160),
                          secs = secs, stringsAsFactors = FALSE)
      }
      if (verbose) cat(sprintf("%-20s %-28s %-9s %5.1fs %s\n", gn, cf$label, row$outcome, secs, row$detail))
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}
