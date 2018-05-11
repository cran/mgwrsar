#' bias_gwrsar
#' to be documented
#' @usage bias_gwrsar(model, data, Model, DGP)
#' @param model to be documented
#' @param data to be documented
#' @param Model to be documented
#' @param DGP to be documented
#' @keywords internal
#' @return to be documented
bias_gwrsar <-
function(model, data, Model, DGP){
    n <- nrow(data)
    names_X = attr(terms(model.frame(model$formula, data)), "term.labels")
    names_X = c("Intercept", names_X)
    data$Intercept = rep(1, n)
    k = length(names_X)
    BETA = matrix(0, ncol = k, nrow = n)
    for (i in 0:(k - 1)) {
        BETA[, i + 1] <- data[, paste("Beta", i, sep = "")]
    }
    BETAV = model$Betav
    BETAC = model$Betac
    kv = ncol(BETAV)
    kc = length(BETAC)
    if (!is.null(model$fixed_vars)) {
        index_XC = which(names_X %in% model$fixed_vars)
        index_XV = (1:k)[-index_XC]
    }
    else {
        index_XC = NULL
        index_XV = (1:k)
    }
    if (Model == "SAR")
        index_XC = 1:(kc - 1)
    if (Model == "OLS")
        index_XC = 1:kc
    if (Model == "MGWRSAR_1_kc_0") {
        index_XC = 1:kc
        index_XV = NULL
    }
    cat("\n######### DGP = ", DGP, " #########\n")
    cat("\n######### model_DGP = ", Model, " #########\n")
    if (Model %in% c("OLS", "GWR", "MGWR"))
        BIAS = matrix(0, ncol = k, nrow = 1)
    else BIAS = matrix(0, ncol = k + 1, nrow = 1)
    if (Model == "MGWRSAR_1_0_kv")
        kv = kv - 1
    if (Model %in% c("MGWRSAR_0_0_kv", "GWR", "MGWRSAR_1_0_kv")) {
        for (i in 1:kv) BIAS[i] <- mean(BETA[, i] - BETAV[, i])
        if (Model == "MGWRSAR_0_0_kv")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betac)
        if (Model == "MGWRSAR_1_0_kv")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betav[length(model$Betav)])
    }
    if (Model %in% c("MGWR", "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv")) {
        for (i in index_XV) {
            j = which(i == index_XV)
            BIAS[i] <- mean(BETA[, i] - BETAV[, j])
        }
        for (i in index_XC) {
            j = which(i == index_XC)
            BIAS[i] <- mean(BETA[, i] - BETAC[j])
        }
        if (Model == "MGWRSAR_0_kc_kv")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betac[length(model$Betac)])
        if (Model == "MGWRSAR_1_kc_kv")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betav[length(model$Betav)])
    }
    if (Model %in% c("MGWRSAR_1_kc_0", "OLS", "SAR")) {
        for (i in index_XC) BIAS[i] <- mean(BETA[, i] - BETAC[i])
        if (Model == "MGWRSAR_1_kc_0")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betav[length(model$Betav)])
        if (Model == "SAR")
            BIAS[ncol(BIAS)] <- mean(data$lambda - model$Betac[length(model$Betac)])
    }
    BIAS
}
