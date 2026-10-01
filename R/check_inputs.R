# INTERNAL: explicit validation of MGWRSAR() inputs
#
# Bad inputs used to surface deep in the code as index, dimension or NA errors
# (campaign of 2026-10). Each check stops with a message naming the argument.

.mgwrsar_models <- c("OLS", "SAR", "GWR", "GWR_glm", "GWR_glmboost",
                     "GWR_gamboost_linearized", "multiscale_gwr", "GWR_multiscale",
                     "MGWR", "MGWRSAR_0_0_kv", "MGWRSAR_0_kc_kv", "MGWRSAR_1_0_kv",
                     "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0")
.mgwrsar_kernels <- c("gauss", "bisq", "epane", "triangle", "tcub", "rectangle",
                      "shepard", "edk", "kdist")

.mgwrsar_check_inputs <- function(data_cols, coords, kernels, H, Model, control, n) {
  if (!(Model %in% .mgwrsar_models))
    stop(sprintf("unknown Model '%s'; use one of %s.", Model,
                 paste(.mgwrsar_models, collapse = ", ")), call. = FALSE)
  Type <- control$Type
  if (!(Type %in% c("GD", "GDT", "T")))
    stop(sprintf("unknown control$Type '%s'; use 'GD' (space), 'GDT' (space-time) or 'T' (time).", Type), call. = FALSE)
  if (anyNA(data_cols))
    stop("the data contain missing values in the variables of the formula: remove or impute them first.", call. = FALSE)
  if (Model != "OLS") {
    if (is.null(coords)) stop("coords must be provided", call. = FALSE)
    if (anyNA(coords)) stop("coords contain missing values.", call. = FALSE)
    if (Type %in% c("GDT", "T")) {
      if (is.null(control$Z))
        stop(sprintf("control$Z (the time variable) is required for Type = '%s'.", Type), call. = FALSE)
      if (anyNA(control$Z)) stop("control$Z (time) contains missing values.", call. = FALSE)
      if (length(control$Z) != n)
        stop(sprintf("control$Z has length %d but the data have %d rows.", length(control$Z), n), call. = FALSE)
    }
    if (!(Model %in% c("SAR")) && !isTRUE(control$searchB)) {
      base <- sub("_.*$", "", kernels)
      bad <- base[!(base %in% .mgwrsar_kernels)]
      if (length(bad))
        stop(sprintf("unknown kernel(s) %s; use one of %s (optionally suffixed by '_past' or a cyclic period).",
                     paste(sQuote(bad), collapse = ", "), paste(.mgwrsar_kernels, collapse = ", ")), call. = FALSE)
      if (!is.null(H)) {
        if (!is.numeric(H) || anyNA(H))
          stop("H (the bandwidths) must be numeric without missing values.", call. = FALSE)
        adaptive <- rep_len(as.logical(if (is.null(control$adaptive)) FALSE else control$adaptive), length(H))
        if (any(!adaptive & !(H > 0)))
          stop(sprintf("a fixed bandwidth must be positive (H = %s).", paste(H[!adaptive & !(H > 0)], collapse = ", ")), call. = FALSE)
        if (any(adaptive & !(H >= 1)))
          stop(sprintf("an adaptive bandwidth is a number of neighbours, at least 1 (H = %s).", paste(H[adaptive & !(H >= 1)], collapse = ", ")), call. = FALSE)
      }
    }
  }
  # the W slot of a fitted model is an empty S4 prototype when the model has no
  # spatial lag; predict() and the bootstrap pass it back as control$W
  W <- control$W
  if (!is.null(W) && (is.matrix(W) || is.data.frame(W) || inherits(W, "Matrix"))) {
    dW <- dim(W)
    if (length(dW) != 2L || dW[1] != n || dW[2] != n)
      stop(sprintf("control$W must be an n x n matrix (n = %d); it is %s.", n,
                   paste(dW, collapse = " x ")), call. = FALSE)
  }
  invisible(TRUE)
}
