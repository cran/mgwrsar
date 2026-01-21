#' Simulate Data Generating Processes (DGP) for Multiscale GWR
#'
#' @description
#' The \code{simu_multiscale} function generates synthetic datasets with spatially varying coefficients
#' based on Data Generating Processes (DGP) proposed in the literature. It supports
#' formulations from Fotheringham et al. (2017), Gao et al. (2021), and Geniaux (2024).
#' It also allows for the introduction of spatial autocorrelation (SAR) in the response variable.
#'
#' @param n Integer. The number of observations to simulate. Default is 1000.
#' @param myseed Integer. Random seed for reproducibility. Default is 1.
#' @param type Character. The type of DGP to use. Options are \code{'FT2017'}, \code{'Gao2021'}, or \code{'GG2024'} (default).
#' @param constant Vector. Indices or names of explanatory variables that should have constant (spatially stationary) coefficients. Default is \code{NULL} (all coefficients vary).
#' @param nuls Vector. Indices of explanatory variables that should have a null effect (coefficient = 0). Default is \code{NULL}.
#' @param lambda Numeric. The spatial autoregressive parameter (rho) if a SAR process is desired. Default is \code{NULL} (no spatial autocorrelation).
#' @param NN Integer. The number of nearest neighbors used to construct the spatial weight matrix \code{W} if \code{lambda} is provided. Default is 4.
#' @param config_beta Character. Configuration of the spatial pattern for Beta coefficients (e.g., \code{'default'}).
#' @param config_snr Numeric. The desired Signal-to-Noise Ratio (SNR). Default is 0.7.
#' @param config_eps Character. The distribution of the error term. Options are \code{'normal'} (default), \code{'unif'}, or \code{'Chi2'}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{mydata}{A data frame with the simulated response variable \code{y} and explanatory variables.}
#'   \item{coords}{A matrix of spatial coordinates.}
#'   \item{Betav}{A matrix of the true spatially varying coefficients.}
#'   \item{W}{The spatial weight matrix (if \code{lambda} is not NULL).}
#' }
#'
#' @references
#' Fotheringham, A. S., Yang, W., & Kang, W. (2017). Multiscale geographically weighted regression (MGWR). \emph{Annals of the American Association of Geographers}, 107(6), 1247-1265.
#'
#' Gao, S., Mei, C. L., Xu, Q. X., & Wang, N. (2021). Non-iterative multiscale estimation for spatial autoregressive geographically weighted regression models. \emph{Entropy}, 25(2), 320.
#'
#' Geniaux, G. (2025). Top-Down Scale Approaches for Multiscale GWR with Locally Adaptive Bandwidths. \emph{Springer Nature}, 2021 LATEX template.
#'
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  library(ggplot2)
#'  library(gridExtra)
#'  library(grid)
#'
#'  # Simulate data using Geniaux (2024) DGP
#'  simu <- simu_multiscale(n = 1000, type = 'GG2024', config_snr = 0.7)
#'  mydata <- simu$mydata
#'  coords <- simu$coords
#'
#'  # Visualizing the spatial patterns of the coefficients
#'  p1 <- ggplot(mydata, aes(x, y, col = Beta1)) + geom_point() + scale_color_viridis_c()
#'  p2 <- ggplot(mydata, aes(x, y, col = Beta2)) + geom_point() + scale_color_viridis_c()
#'  p3 <- ggplot(mydata, aes(x, y, col = Beta3)) + geom_point() + scale_color_viridis_c()
#'  p4 <- ggplot(mydata, aes(x, y, col = Beta4)) + geom_point() + scale_color_viridis_c()
#'
#'  grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2,
#'               top = textGrob("DGP Geniaux (2024)", gp = gpar(fontsize = 20, font = 3)))
#' }
#' @export
simu_multiscale <- function(n = 1000, myseed = 1, type = 'GG2024', constant = NULL, nuls = NULL, lambda = NULL, NN = 4, config_beta = 'default', config_snr = 0.7, config_eps = 'normal') {

  # 1. RNG Initialization (Fixed Seed)
  set.seed(myseed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")

  W <- NULL

  # --- DETERMINISTIC SNR FUNCTION ---
  get_snr_fixed <- function(k, config_snr, XB, eps_raw) {
    eps_scaled <- k * eps_raw
    # Inverse SNR formula: finding k such that signal/total ratio matches config_snr
    # (Note: ensure this formula is indeed the one you intend to keep)
    abs(config_snr - (sum(XB^2)) / sum((XB + eps_scaled)^2))
  }

  # 2. Model Type Configuration
  if (config_beta == 'All_B3') type = 'All_B3'
  if (config_beta == 'B4_nonconstant') constant = c(1, 2, 3)
  if (config_beta == 'B1_constant') constant = 1
  if (config_beta == 'B2_constant') constant = 2
  if (config_beta == 'B3_constant') constant = 3
  if (config_beta == 'B4_constant') constant = 4
  if (config_beta == 'B1_nul') nuls = 1
  if (config_beta == 'B2_nul') nuls = 2
  if (config_beta == 'B3_nul') nuls = 3
  if (config_beta == 'B4_nul') nuls = 4
  if (config_beta == 'GG2024_univ') type = 'GG2024_univ'
  if (config_beta == 'low_h') type = 'GG2024_low'
  if (config_beta == 'high_h') type = 'GG2024_high'
  if (config_beta == 'FT2017') type = 'FT2017'
  if (config_beta == 'Gao2021') type = 'Gao2021'
  if (config_beta == 'GG2024_2') type = 'GG2024_2'
  if (config_beta == 'GG2024_8') type = 'GG2024_8'
  if (config_beta == 'GG2024_16') type = 'GG2024_16'
  if (config_beta == 'GG2024_32') type = 'GG2024_32'
  if (config_beta == 'spatiotemp') type = 'spatiotemp'
  if (config_beta == 'spatiotemp_old') type = 'spatiotemp_old'
  if (config_beta == 'spatiotemp_fake') type = 'spatiotemp_fake'
  if (config_beta == 'same_h') type = 'same_h'
  if (config_beta == 'same_ht') type = 'same_ht'

  # 3. Coordinates and Base Variables Generation
  x <- runif(n)
  y <- runif(n)
  coords <- cbind(x, y)

  if (!is.null(lambda)) {
    W <- kernel_matW(H = NN, kernels = 'rectangle', coords = coords, NN = NN, adaptive = TRUE, diagnull = TRUE)
    if (length(lambda) == 1) lambda <- rep(lambda, n)
    iW <- solve(Diagonal(n, 1) - lambda * W)
  } else iW <- NULL

  X1 <- runif(n, -sqrt(3), sqrt(3))
  X2 <- rnorm(n)
  X3 <- rnorm(n) # Default (will be overwritten or used depending on case)

  # --- RAW NOISE GENERATION (FROZEN NOISE) ---
  # Done HERE, once, to ensure Mac/Linux identity
  if (config_eps == 'normal') {
    eps_raw <- rnorm(n, 0, 1)
  } else if (config_eps == 'unif') {
    eps_raw <- runif(n, min = -3, max = 3)
  } else if (config_eps == 'Chi2') {
    eps_raw <- 0.5 * rchisq(n, df = 2, ncp = 0) - 1
  }

  # 4. Beta Calculation According to Type
  Beta1 <- 3 * (x + y) - 1

  if (type == 'FT2017') {
    Beta1 <- 3
    Beta2 <- 1 + (x * 25 + y * 25) / 12
    Beta3 <- 1 + 1 / 324 * (36 - (6 - x * 25 / 2)^2) * (36 - (6 - y * 25 / 2)^2)

  } else if (type == 'Gao2021') {
    Beta2 <- 4 * sin(sqrt(12 * (x - 0.5)^2 + 12 * (y - 0.5)^2))
    Beta3 <- 84 * (x * y) * (1 - x) * (1 - y)

  } else if (type == 'GG2024_univ') {
    Beta1 <- 3 * x - 1
    Beta2 <- 84 * (x / 2) * (1 - x) * 0.5
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - 0.25)^2))) + 2 * (x - 0.5)

  } else if (type == 'GG2024_low') {
    Beta1 <- 3 * x - 1
    Beta2 <- 84 * (x * y) * (1 - x) * (1 - y)
    Beta3 <- 3 * (x - y) - 1
    Beta4 <- 3 * (y - x) - 1

  } else if (type == 'GG2024_high') {
    Beta1 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x + y)^4)))
    Beta2 <- 84 * (x * y) * (1 - x) * (1 - y)
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - y)^2))) + 2 * (x - 0.5)

  } else if (type == 'All_B3') {
    Beta1 <- 3
    Beta2 <- Beta4 <- Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))

  } else if (type == 'spatiotemp') {
    time <- sample(1:365, n, replace = TRUE)
    year <- sample(1:4, n, replace = TRUE)
    x_cor <- x + abs(time - 210) / 210
    x_cor <- x_cor / max(x_cor)
    y_cor <- y + abs(time - 210) / 210
    y_cor <- (y_cor - min(y_cor)) / max(y_cor - min(y_cor))

    direction <- abs((time - 182.5) / 182.5)
    Beta1 <- (1 - direction) * (3 * (x + y) - 1) + direction * (3 * (-x + y))
    Beta1 <- Beta1 - min(Beta1)
    Beta2 <- (84 * (x_cor * y_cor) * (1 - x_cor) * (1 - y_cor))
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - y)^2))) + 2 * (x - 0.5)

  } else if (type == "spatiotemp_old") {
    ratiotime <- 1.5
    time <- seq(0, 1, length.out = n)
    timeratio <- 1 + time * ratiotime
    Beta1 <- (3 * (x + y) - 1) * timeratio
    Beta2 <- (84 * (x * y) * (1 - x) * (1 - y)) * timeratio
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - y)^2))) + 2 * (x - 0.5)

  } else if (type == 'spatiotemp_fake') {
    time <- sample(1:365, n, replace = TRUE)
    year <- sample(1:4, n, replace = TRUE)

    Beta1 <- 3 * (x + y) - 1
    Beta2 <- (84 * (x * y) * (1 - x) * (1 - y))
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - y)^2))) + 2 * (x - 0.5)

  } else if (type == 'same_ht') {
    time <- sample(1:365, n, replace = TRUE)
    year <- sample(1:4, n, replace = TRUE)

    x_cor <- x + abs(time - 210) / 210
    x_cor <- x_cor / max(x_cor)
    y_cor <- y + abs(time - 210) / 210
    y_cor <- (y_cor - min(y_cor)) / max(y_cor - min(y_cor))

    z1 <- x_cor - 0.2
    z2 <- y_cor - 0.2
    z3 <- x_cor + 0.2
    z4 <- y_cor + 0.2
    Beta1 <- 84 * (z1 * z2) * (1 - z1) * (1 - z2)
    Beta2 <- 84 * (z3 * z4) * (1 - z3) * (1 - z4)
    Beta3 <- 84 * (z1 * z4) * (1 - z1) * (1 - z4)
    Beta4 <- 84 * (z3 * z2) * (1 - z3) * (1 - z2)

  } else if (type == 'same_h') {
    time <- sample(1:365, n, replace = TRUE)
    year <- sample(1:4, n, replace = TRUE)

    z1 <- x - 0.2
    z2 <- y - 0.2
    z3 <- x + 0.2
    z4 <- y + 0.2
    Beta1 <- 84 * (z1 * z2) * (1 - z1) * (1 - z2)
    Beta2 <- 84 * (z3 * z4) * (1 - z3) * (1 - z4)
    Beta3 <- 84 * (z1 * z4) * (1 - z1) * (1 - z4)
    Beta4 <- 84 * (z3 * z2) * (1 - z3) * (1 - z2)

  } else { # DEFAULT
    Beta2 <- 84 * (x * y) * (1 - x) * (1 - y)
    Beta3 <- 4 * x * sin(sqrt((6^2 * (-0.4 + x * 2)^4)))
    Beta4 <- 4 * sin(sqrt((6^3 * (-x - y)^2))) + 2 * (x - 0.5)
  }

  # Constant Coefficients Handling
  if (length(constant) == 1) {
    if (constant == 1) Beta1 = 3
    if (constant == 2) Beta2 = 3
    if (constant == 3) Beta3 = 3
    if (constant == 4) Beta4 = 3
  } else {
    for (i in 1:length(constant)) assign(paste0('Beta', constant[i]), 3)
  }

  # Null Coefficients Handling
  if (length(nuls) == 1 & length(constant) == 0) {
    if (nuls == 1) Beta1 = 0
    if (nuls == 2) Beta2 = 0
    if (nuls == 3) Beta3 = 0
    if (nuls == 4) Beta4 = 0
  } else {
    for (i in 1:length(nuls)) assign(paste0('Beta', nuls[i]), 0)
  }

  # 5. Construction of Y and mydata (Robustification)
  if (type %in% c('FT2017', 'Gao2021')) {
    XB <- Beta1 + Beta2 * X1 + Beta3 * X2

    # SNR Optimization with frozen noise
    k <- (optimize(get_snr_fixed, lower = 1, upper = 100, config_snr = config_snr, XB = XB, eps_raw = eps_raw))$minimum
    eps <- k * eps_raw

    Y <- XB + eps
    if (!is.null(iW)) Y <- as.numeric(iW %*% Y)
    mydata <- data.frame(Y, X1, X2, Beta1, Beta2, Beta3, eps, x, y)

  } else if (type %in% c('GG2024_univ')) {
    # For this specific case, you had a different noise logic
    if (config_snr >= 0.9) eps <- rnorm(n, 0, 0.1)
    else if (config_snr > 0.7) eps <- rnorm(n, 0, 0.25)
    else if (config_snr > 0.5) eps <- rnorm(n, 0, 0.5)
    else eps <- eps_raw # Fallback

    Y <- Beta1 + Beta2 * X1 + Beta3 * X2 + Beta4 * X3 + eps
    Y1 <- Beta1 + eps
    Y2 <- Beta2 * X1 + eps
    Y3 <- Beta3 * X2 + eps
    Y4 <- Beta4 * X3 + eps
    if (!is.null(iW)) Y <- as.numeric(iW %*% Y)
    mydata <- data.frame(Y, Y1, Y2, Y3, Y4, X1, X2, X3, Beta1, Beta2, Beta3, Beta4, eps, x, y)

  } else if (type == 'GG2024_2') {
    XB <- Beta1 + Beta2 * X1
    k <- (optimize(get_snr_fixed, lower = 1, upper = 100, config_snr = config_snr, XB = XB, eps_raw = eps_raw))$minimum
    eps <- k * eps_raw
    Y <- XB + eps
    mydata <- data.frame(Y, X1, Beta1, Beta2, eps, x, y)

  } else if (type == 'GG2024_8') {
    X4 <- runif(n, -sqrt(3), sqrt(3))
    X5 <- rnorm(n)
    X6 <- runif(n, -sqrt(3), sqrt(3))
    X7 <- rnorm(n)
    Beta5 <- 1 + (x * 25 + y * 25) / 12
    Beta6 <- 1 + 1 / 324 * (36 - (6 - x * 25 / 2)^2) * (36 - (6 - y * 25 / 2)^2)
    Beta7 <- 4 * sin(sqrt(12 * (x - 0.5)^2 + 12 * (y - 0.5)^2))
    Beta8 <- 84 * (x * y) * (1 - x) * (1 - y)
    XB <- Beta1 + Beta2 * X1 + Beta3 * X2 + Beta4 * X3 + Beta5 * X4 + Beta6 * X5 + Beta7 * X6 + Beta8 * X7
    k <- (optimize(get_snr_fixed, lower = 1, upper = 100, config_snr = config_snr, XB = XB, eps_raw = eps_raw))$minimum
    eps <- k * eps_raw
    Y <- XB + eps
    if (!is.null(iW)) Y <- as.numeric(iW %*% Y)
    mydata <- data.frame(Y, X1, X2, X3, X4, X5, X6, X7, Beta1, Beta2, Beta3, Beta4, Beta5, Beta6, Beta7, Beta8, eps, x, y)

  } else if (type == 'GG2024_16') {
    # (Creation of X4 to X15...)
    X4 <- runif(n, -sqrt(3), sqrt(3)); X5 <- rnorm(n); X6 <- runif(n, -sqrt(3), sqrt(3)); X7 <- rnorm(n)
    X8 <- runif(n, -sqrt(3), sqrt(3)); X9 <- rnorm(n); X10 <- runif(n, -sqrt(3), sqrt(3)); X11 <- rnorm(n)
    X12 <- runif(n, -sqrt(3), sqrt(3)); X13 <- rnorm(n); X14 <- runif(n, -sqrt(3), sqrt(3)); X15 <- rnorm(n)

    Beta5 <- 1 + (x * 25 + y * 25) / 12
    Beta6 <- 1 + 1 / 324 * (36 - (6 - x * 25 / 2)^2) * (36 - (6 - y * 25 / 2)^2)
    Beta7 <- 4 * sin(sqrt(12 * (x - 0.5)^2 + 12 * (y - 0.5)^2))
    Beta8 <- 84 * (x * y) * (1 - x) * (1 - y)
    Beta9 <- Beta1; Beta10 <- Beta2; Beta11 <- Beta3; Beta12 <- Beta4
    Beta12 <- Beta4 <- 3 # Remove as in your original code
    Beta13 <- Beta5; Beta14 <- Beta6; Beta15 <- Beta7; Beta16 <- Beta8

    XB <- Beta1 + Beta2 * X1 + Beta3 * X2 + Beta4 * X3 + Beta5 * X4 + Beta6 * X5 + Beta7 * X6 + Beta8 * X7 +
      Beta9 * X8 + Beta10 * X9 + Beta11 * X10 + Beta12 * X11 + Beta13 * X12 + Beta14 * X13 + Beta15 * X14 + Beta16 * X15

    k <- (optimize(get_snr_fixed, lower = 1, upper = 100, config_snr = config_snr, XB = XB, eps_raw = eps_raw))$minimum
    eps <- k * eps_raw
    Y <- XB + eps
    if (!is.null(iW)) Y <- as.numeric(iW %*% Y)
    mydata <- data.frame(Y, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, Beta1, Beta2, Beta3, Beta4, Beta5, Beta6, Beta7, Beta8, Beta9, Beta10, Beta11, Beta12, Beta13, Beta14, Beta15, Beta16, eps, x, y)

  } else {
    # Default case
    XB <- Beta1 + Beta2 * X1 + Beta3 * X2 + Beta4 * X3
    k <- (optimize(get_snr_fixed, lower = 1, upper = 100, config_snr = config_snr, XB = XB, eps_raw = eps_raw))$minimum
    eps <- k * eps_raw
    Y <- XB + eps
    if (!is.null(iW)) Y <- as.numeric(iW %*% Y)
    mydata <- data.frame(Y, X1, X2, X3, Beta1, Beta2, Beta3, Beta4, eps, x, y)
  }

  mydata$Intercept <- 1
  if (config_beta %in% c('spatiotemp', 'spatiotemp_old', 'spatiotemp_fake', 'same_h', 'same_ht')) {
    mydata$time <- time
    if (config_beta != 'spatiotemp_old') mydata$year <- year
  }
  if (!is.null(lambda)) mydata$lambda <- lambda

  list(mydata = mydata, coords = coords, W = W)
}

