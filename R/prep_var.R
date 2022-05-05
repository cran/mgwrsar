#' prep_var
#' to be documented
#' @usage prep_var(gwrenv)
#' @param gwrenv to be documented
#' @noRd
#' @return to be documented
prep_var<-function(gwrenv){
  if (is.null(gwrenv$coord)) {
    if (class(gwrenv$data) %in% c("SpatialPointsDataFrame", "SpatialGridDataFrame","SpatialPixelsDataFrame")) gwrenv$coord = as.matrix(coordinates(gwrenv$data)) else stop("coord required")
  }

  if (length(gwrenv$kernels) > 1) gwrenv$S = as.matrix(cbind(gwrenv$coord, gwrenv$Z)) else gwrenv$S = as.matrix(gwrenv$coord)
  if(gwrenv$Model!='SAR'){
  if (!gwrenv$searchB  & (is.null(gwrenv$H[1]) | is.null(gwrenv$kernels[1]))) stop("kernels list and bandwidths H required")
  if (!gwrenv$searchB & gwrenv$adaptive[1] & gwrenv$kernels[1]!='gauss') {
    gwrenv$NN=gwrenv$H[1]+2
  }
  if (is.null(gwrenv$fixed_vars) & gwrenv$Model %in% c("MGWR", "MGWRSAR_0_kc_kv","MGWRSAR_1_kc_kv"))
    stop("You must provide fixed_vars for mixed models")
}
  if (is.null(gwrenv$W) & !gwrenv$searchB & gwrenv$Model %in% c("SAR", "MGWRSAR_1_0_kv", "MGWRSAR_0_0_kv","MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0"))
    stop("You must provide W for models with spatial dependence")
  if (!is.null(gwrenv$fixed_vars) & gwrenv$Model %in% c("GWR", "SAR", "MGWRSAR_1_0_kv","MGWRSAR_0_0_kv")) {
    gwrenv$fixed_vars = NULL
    if (gwrenv$verbose)
      cat("\n-----------------------------------------------------\nfixed_vars set to NULL because model= ",
          gwrenv$Model, "\n-----------------------------------------------------\n")
  }
  if (!is.null(gwrenv$W) & gwrenv$Model %in% c("GWR", "OLS", "MGWR")) {
    if (gwrenv$verbose)
      cat("\n-----------------------------------------------------\nW not used because model= ",
          gwrenv$Model, "\n-----------------------------------------------------\n")
  }
  gwrenv$mf <- model.frame(gwrenv$formula, gwrenv$data)
  if(!is.null(gwrenv$new_data)) gwrenv$new_mf <- model.frame(gwrenv$formula, gwrenv$new_data)
  gwrenv$mt <- attr(x = gwrenv$mf, which = "terms")
  gwrenv$X = model.matrix(object = gwrenv$mt, data = gwrenv$mf, contrasts.arg = gwrenv$contrasts)
  if(!is.null(gwrenv$new_data)) gwrenv$new_X = model.matrix(object = gwrenv$mt, data = gwrenv$new_mf, contrasts.arg = gwrenv$contrasts)
  gwrenv$Y <- model.extract(gwrenv$mf, "response")
  idx1 <- match("(Intercept)", colnames(gwrenv$X))
  if(!is.null(gwrenv$new_data)) new_idx1 <- match("(Intercept)", colnames(gwrenv$new_X))

  if (!is.na(idx1))
    colnames(gwrenv$X)[idx1] <- "Intercept"
  if(!is.null(gwrenv$new_data)) {
    if (!is.na(new_idx1))
      colnames(gwrenv$new_X)[idx1] <- "Intercept"
  }
  if (!is.null(gwrenv$fixed_vars)) {
    idx.fixed <- as.numeric(na.omit(match(gwrenv$fixed_vars, colnames(gwrenv$X))))
    gwrenv$XC <- as.matrix(gwrenv$X[, idx.fixed])
    colnames(gwrenv$XC) <- colnames(gwrenv$X)[idx.fixed]
    if (length(idx.fixed) < ncol(gwrenv$X)) {
      gwrenv$XV <- as.matrix(gwrenv$X[, -idx.fixed])
      colnames(gwrenv$XV) <- colnames(gwrenv$X)[-idx.fixed]
    } else gwrenv$XV = NULL
    if(!is.null(gwrenv$new_data)) {
      gwrenv$new_XC <- as.matrix(gwrenv$new_X[, idx.fixed])
      colnames(gwrenv$new_XC) <-colnames(gwrenv$new_X)[idx.fixed]
      if (length(idx.fixed) < ncol(gwrenv$X)) {
        gwrenv$new_XV <- as.matrix(gwrenv$new_X[, -idx.fixed])
        colnames(gwrenv$new_XV) <- colnames(gwrenv$new_X)[-idx.fixed]
      } else gwrenv$new_XV = NULL
    }
  }
  else {
    gwrenv$XV = as.matrix(gwrenv$X)
    gwrenv$XC = NULL
    if(!is.null(gwrenv$new_data)) {
      gwrenv$new_XV = as.matrix(gwrenv$new_X)
      gwrenv$new_XC = NULL
    }
  }
  gwrenv$coord = as.matrix(gwrenv$coord)
  # if (is.null(gwrenv$W))
  #   gwrenv$W <- as(Matrix(0, nrow = gwrenv$n, ncol = gwrenv$n), "dgCMatrix")
  gwrenv$names_betac = colnames(gwrenv$XC)
  gwrenv$names_betav = colnames(gwrenv$XV)
  if (gwrenv$Model %in% c("OLS"))
    gwrenv$names_betac = colnames(gwrenv$X)
  if (gwrenv$Model %in% c("SAR"))
    gwrenv$names_betac = c(colnames(gwrenv$X), "lambda")
  if (gwrenv$Model %in% c("MGWRSAR_0_kc_kv", "MGWRSAR_0_0_kv"))
    gwrenv$names_betac = c(gwrenv$names_betac, "lambda")
  if (gwrenv$Model %in% c("MGWRSAR_1_0_kv", "MGWRSAR_1_kc_kv"))
    gwrenv$names_betav = c(gwrenv$names_betav, "lambda")
  if (gwrenv$Model == "MGWRSAR_1_kc_0") {
    gwrenv$names_betav = c("lambda")
    gwrenv$names_betac = colnames(gwrenv$X)
  }
  gwrenv$MykernelS = gwrenv$kernels
  gwrenv$HH = gwrenv$H
  gwrenv$Y = as.matrix(gwrenv$Y)
  gwrenv$X = as.matrix(gwrenv$X)
  gwrenv$ALL_X = as.matrix(gwrenv$X)
  if (!is.null(gwrenv$XC))
    gwrenv$XC = as.matrix(gwrenv$XC)
  if (!is.null(gwrenv$XV))
    gwrenv$XV = as.matrix(gwrenv$XV)
  if (is.null(gwrenv$TP)) gwrenv$TP=1:gwrenv$n
  if (!is.null(gwrenv$S_out)) gwrenv$TP=1:nrow(gwrenv$S_out)
  gwrenv
}
