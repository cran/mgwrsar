NEWS/ChangeLog
-----------------------------
# mgwrsar 1.4.4 (2026-10-28)

Released on CRAN.

## Bug fixes

- Local standard errors (`@sev`) of `MGWRSAR()` (GWR, mixed GWR, SAR variants) and of `TDS_MGWR()` without `get_AIC` now use the variance of the local smoother with a local residual variance. They ignored the effective local sample size of the normalised kernel weights, and the fixed part of the residual in mixed models, so that they were too small (2.8 to 4.5 times on a 2,667-observation `tds_mgtwr` fit). t-values, p-values and significance maps change. Affected: all previous versions. Coefficients, fitted values, hat matrices, AICc and the standard errors of `TDS_MGWR(get_AIC = TRUE)` are unchanged.
- `predict(method_pred = "TP")` now uses the control of the training call (kernel type, adaptive bandwidth, `W`), and for a mixed model (`Model = "MGWR"`) keeps the fixed coefficients and re-estimates the varying ones on the partial residual. Affected: 1.3.1 to 1.4.3, where the defaults were used instead, and every earlier version for mixed models.
- New tests compare every prediction mode and the local standard errors with references computed by hand.

# mgwrsar 1.4.3 (2026-10-27)

Not released on CRAN.

## Changed behaviour

- `TDS_MGWR()` returns the terminal sweep of the backfitting, as MGWR backfitting does, instead of the sweep with the smallest in-sample RMSE (1.3.1 to 1.4.2). `control_tds$retained_state = "best"` restores the former behaviour.
- `TDS_MGWR()`: after the search, the coefficients are refitted at fixed bandwidths until convergence (`control_tds$polish`, default 3 sweeps; 0 disables it).
- `TDS_MGWR()`: the safeguard against saturated candidate fits is a unitless ratio of residual variances (`control_tds$saturation_ratio`, default `exp(-2)`), with counts in the new slot `saturation`. Up to 1.4.2 it depended on the unit of the response and on `get_AIC`.
- `TDS_MGWR()`: the importance-driven update order uses the standard deviation of `beta_k x_k` and no longer depends on the units or the level of the covariates (up to 1.4.2 it did). The path of the descent may change.
- `predict()` for `tds_mgwr` / `tds_mgtwr`: a target with very small prediction weights no longer gets a non-finite prediction, and a target without any positive weight takes the coefficient of its nearest training point, with a warning.
- New levers of `control_tds`: `ordering = "random"`, besides `"cyclic"` and the default `"Importance"`; `importance_order` (default `"decreasing"`); `update_last` (default `NULL`, no covariate forced last).

# mgwrsar 1.4.2 (2026-10-06)

Not released on CRAN.

## Bug fixes

- `TDS_MGWR()` with pinned bandwidths, an adaptive kernel and `NN < n`: a coefficient pinned at the global value is fitted as a global component, as `predict()` assumes. Affected: 1.4.1.
- `TDS_MGWR(control = list(SE = TRUE))`: `@sev` and `@edf_k` are those of the returned sweep, not of the last one. Affected: 1.4.1 (relative errors up to about 18%).
- `search_bandwidths()` no longer discards the true minimum of a grid round, which could return a worse local minimum. Affected: 1.4.1. The AICc of a saturated binomial GWR is now `Inf`.

# mgwrsar 1.4.1 (2026-10-01)

Released on CRAN.

## New features

- New `gtwr_HWB2010()`: the GTWR of Huang, Wu and Barry (2010) with its `(h_ST, tau)` parameterisation, and converters `bw_hwb2gdt()`, `bw_gdt2hwb()` and `as_gtwr()`.
- `summary()` gains `fdr_method = "spatial_BY"`, a Benjamini-Yekutieli correction calibrated on the bandwidth.
- New `control_tds` levers: per-coefficient and per-axis pinned bandwidths (`H`, `Ht`, now working for `Type = 'GDT'`), a grid in distances (`V_dist`), a floor of the spatial grid for panel data (`panel_floor`) and `TRUEBETA` for simulations.

## Performance

- Native engines for GWR, mixed GWR and MGWR, native kernel weights, a session cache of weight matrices and truncated compact kernels: `TDS_MGWR()` is about 4 to 5 times faster than in 1.3.2 and `search_bandwidths()` about 1.5 to 2 times, with unchanged results.
- SMUT is no longer imported; dense products use Eigen or the BLAS of R (`options(mgwrsar.matprod = )`).

## Bug fixes and robustness

- The AICc is `Inf` for a saturated fit, which previous versions rewarded with a penalty of the wrong sign.
- Local systems are solved by a column-pivoted QR: nearly singular designs, constant or collinear predictors and isolated points no longer give non-finite coefficients or leverages outside [0, 1].
- `TDS_MGWR()`: after a sweep that raises the RMSE the descent continues from that sweep; 1.3.2 restarted from an inconsistent state that froze the descent.
- The random number generator of the session is restored on exit by `MGWRSAR()`, `search_bandwidths()`, `TDS_MGWR()` and the other functions that draw internally. Affected: previous versions.
- `search_bandwidths()` no longer closes every connection of the session, and its golden-section refinement can no longer loop forever.
- Explicit errors replace index or dimension errors for invalid inputs (unknown model, kernel or criterion, missing values, invalid bandwidths, `W` or `fixed_vars`).
- `fitted()` is exported for `mgwrsar` objects; vignettes and `citation()` render correctly; new adversarial test campaign.

# mgwrsar 1.3.2 (2026-03-02)

Released on CRAN.

## TDS algorithms

- Introducing `tds_mgtwr` model for Multiscale GTWR with additive or multiplicative (cyclic or acyclic) spatio-temporal kernels.
- Importance-driven update schedule that prioritizes covariates according to their current scale-normalized contribution to the fitted signal.
- Fixed several edge cases in `TDS_MGWR()`.

## Control-parameter safety

- Added stronger guards on neighborhood and search controls.
- `NN` is now capped to `n` in both `MGWRSAR()` and `TDS_MGWR()` to avoid invalid k-NN requests.
- In `TDS_MGWR()`, `control_tds$nns` is now capped to `min(round(n/8), round(maxit/2))`; if a user-provided value is too large, it is truncated with a warning.

## Parallel search reliability

- Improved stability of parallel bandwidth-search workflows (foreach-based execution), with safer fallback behavior.

## Testing and maintenance

- Added/extended tests for MGTWR/TDS-related scenarios.

# mgwrsar 1.3.1 (2026-01-21)

Released on CRAN.

## Major performance improvements

- Complete rewrite of the local fitting engine using RcppArmadillo (pivotal QR + banded optimizations).
  This substantially accelerates GWR, mixed-GWR and MGWR estimations, with average speed-ups of ×4 and memory usage reduced by 30–50%.

## New bandwidth selection engine

- New function search_bandwidth(), a unified wrapper for 1D (space) and 2D (space and time) bandwidth optimisation.
  It supports multi-round grid refinement, integrates golden-section search, and relies on forking for parallel evaluation.

## Visualisation

- New interactive plot methods based on plotly, allowing dynamic exploration of local coefficients, bandwidth paths, and model diagnostics.

## TDS algorithms

- Significant improvements to tds_mgwr and tds_mgtwr models:
- smoother and more stable AICc-based decisions during bandwidth boosting,
- refined sequential optimisation for spatial and spatio-temporal kernels,
- improved handling of edge cases and isolated observations,
- new internal diagnostics for convergence monitoring.

## Improved handling of isolated points

- Better detection and fallback to OLS for observations receiving zero weight under non-adaptive kernels—avoiding silent numerical instabilities.

## Parallelisation robustness

- More reliable parallel execution:
- cleanup on interruptions,
- fallback to sequential execution when requested.

## Predictive methods

- More stable predict_mgwrsar() logic with safer handling of model@mycall$control, avoiding previous errors when called immediately after bandwidth optimisation.

## Improved numerical stability

- Several fixes related to:
- QR pivoting in local regressions,
- normalization of kernel weights,
- avoidance of underflow in Gaussian kernels for very small bandwidths.

## Cross-platform build stability

- Fixes ensuring compatibility on macOS ARM, Linux (GCC ≥12), and Windows Rtools; better BLAS thread control (OPENBLAS, MKL, VECLIB).

## Reproducibility across platforms

- Adoption of the L’Ecuyer–CMRG random number generator with inversion-based normal deviates, ensuring bitwise-stable stochastic behaviour across all platforms and parallel backends.

# mgwrsar 1.2 (2025-06-24)

Not released on CRAN.

* Introducing tds_mgtwr model allowing to estimate spatio-temporal multiscale GWR (MGTWR)
* Introducing golden_search_2d_bandwidth for automatic bandwidth selection with spatio-temporal GWR (GTWR)

# mgwrsar 1.1 (2024-12-24)

Released on CRAN.
* Introducing top-down scale/multiscale GWR (tds_mgwr), adaptive top-down scale/multiscale GWR (atds_mgwr) and regular multiscale GWR (multiscale_gwr)
* Introducing generalized GWR for binomial (bionomial and quasibinomial families).
* Introducing GAM/GWR with gradient descent boosting.
* Improving kernels for spatio-temporal GWR to introduce seasonnality (GDT)
* Removing unused experimental General kernel Product functions (GDX,GDC)
* Correcting a bug for MGWRSAR_x_x_x models with spatial autocorrelation when SE=TRUE.
* Correcting a bug for GWR model with a single explanatory variable when doMC=TRUE.
* Correcting a bug for formula without explicit names of variable (like "~.").
* Correcting a bug in GWR with parallel computation (split error in gwr_beta)

# mgwrsar 1.0.5 (2023-11-16)

Released on CRAN.
* Removing dependency to qlcMatrix
* Introducing experimental multiscale GWR Model
* Introducing experimental GWR with glm family
* Introducing experimental GWR with glmboost Model
* Adding the ability to calculate the trace of S without calculating the standard deviation matrix (control$get_ts parameter)
* bandwidths_mgwrsar subroutines improved
* AICc criteria added for bandwidth search
* Removing the remove_local_outlier method
* Rename all variables "coord" to "coords"
* Rename the coordinates in the data example to c('x','y')

# mgwrsar 1.0.4 (2023-03-01)

Released on CRAN.
* Improved kernel_matW function by reintroducing possibility of island with no neighbours for non adaptie kernel
* Introduced an experimental version of a backfitting algorithm for the estimation of the Multiscale GWR with a selection of the bandiwth of each variable by cross validation (LOOCV). The predict_mgwrsar function allows to make predictions from a model of the 'multiscale_gwr' class.

# mgwrsar 1.0.3 (2020-11-18)

Not released on CRAN.
* Fixed a bug for plot_mgwrsar
* Fixed a bug for computing SE and edf when system is computationally singular
* Fixed a bug for the function plot_mgwrsar coming from the way of naming the columns of the coordinates matrix returned in the MGWRSAR function.
* Added a warning for cases where isgcv=TRUE and SE=TRUE when calling the MGWRSAR function.
* Improvement of plot_effects function to plot also effect of spatially non-varying parameters
* Improvement of plot_mgwrsar function to allow control of leaflet tile and number of quantile in legend.
* WARNING:  function to predict with target points using TP options have to be checked for model with spatial autocorrelation.

# mgwrsar 1.0.2 (2020-06-01)

Not released on CRAN.

* Fixed a bug in the kernel_matW function when there are duplicted spatal coordinates: duplicated coordinates are jittered with a warning if this is the case.

# mgwrsar 1.0.1 (2020-05-04)

Not released on CRAN.

* All models of mgwrsar package based on local linear regression can now be estimated using a target points set. Several functions that allows to choose an optimal set of target points to obtain a faster approximation of GWR coefficients has been added.

* Predictions on new data can now be done using the jacknife estimation method instead of spatial extrapolation of local coefficients from a preliminarily estimated model. Only the optimal value of the bandwidth modeled from the initial data is then used. In the function 'predict_mgwrsar', if the parameter method_pred = 'TP' (default), the prediction is done by recalculating a MGWRSAR model with the new data as target points keeping the bandwidth at the optimal value chosen with the training data, otherwise if method_pred= ('tWtp_model', 'model', 'shepard') then a matrix is used for the spatial extrapolation of the estimated coefficients, and prediction are done using these extrapolated coefficients (as in the previous version of mgwrsar package).

* 'KNN' function is deprecated and replaced by 'kernel_matW' function that allows to build spatial weight matrix and interaction matrix based on General Kernel Product. In kernel_matW function it's possible to specify the maximum number of neighbors to consider in gaussian kernel (rough gaussian kernel) to increase speed and sparsity of weights matrix.

* Fast computation of local OLS coefficients in the previous version of mgwrsar package (0.1) uses non pivotal computation that may provides undesirable results in presence of strong colinearity. The RCCP 'fastlmLLT_C' function has been replaced by R native lm.fit function in this realease.

* A new ploting function has been added: plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables

# mgwrsar 0.1 (2018-05-11)

Released on CRAN.
* First release.

