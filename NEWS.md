NEWS/ChangeLog
-----------------------------

# 1.3.1 2026-01-21
Major performance improvements
	•	Complete rewrite of the local fitting engine using RcppArmadillo (pivotal QR + banded optimizations).
This substantially accelerates GWR, mixed-GWR and MGWR estimations, with average speed-ups of ×4 and memory usage reduced by 30–50%.

New bandwidth selection engine
	•	New function search_bandwidth(), a unified wrapper for 1D (space) and 2D (space and time) bandwidth optimisation.
It supports multi-round grid refinement, integrates golden-section search, and relies on forking for parallel evaluation.

Visualisation
	•	New interactive plot methods based on plotly, allowing dynamic exploration of local coefficients, bandwidth paths, and model diagnostics.

TDS algorithms
	•	Significant improvements to tds_mgwr and tds_mgtwr models:
	•	smoother and more stable AICc-based decisions during bandwidth boosting,
	•	refined sequential optimisation for spatial and spatio-temporal kernels,
	•	improved handling of edge cases and isolated observations,
	•	new internal diagnostics for convergence monitoring.

Improved handling of isolated points
	•	Better detection and fallback to OLS for observations receiving zero weight under non-adaptive kernels—avoiding silent numerical instabilities.

Parallelisation robustness
	•	More reliable parallel execution:
	•	cleanup on interruptions,
	•	fallback to sequential execution when requested.

Predictive methods
	•	More stable predict_mgwrsar() logic with safer handling of model@mycall$control, avoiding previous errors when called immediately after bandwidth optimisation.

Improved numerical stability
	•	Several fixes related to:
	•	QR pivoting in local regressions,
	•	normalization of kernel weights,
	•	avoidance of underflow in Gaussian kernels for very small bandwidths.

Cross-platform build stability
	•	Fixes ensuring compatibility on macOS ARM, Linux (GCC ≥12), and Windows Rtools; better BLAS thread control (OPENBLAS, MKL, VECLIB).
	
Reproducibility across platforms
	•	Adoption of the L’Ecuyer–CMRG random number generator with inversion-based normal deviates, ensuring bitwise-stable stochastic behaviour across all platforms and parallel backends.

	
# 1.2 2025-6-24 (unreleased version)
* Introducing tds_mgtwr model allowing to estimate spatio-temporal multiscale GWR (MGTWR)
* Introducing golden_search_2d_bandwidth for automatic bandwidth selection with spatio-temporal GWR (GTWR)

# 1.1 2024-12-24
* Introducing top-down scale/multiscale GWR (tds_mgwr), adaptive top-down scale/multiscale GWR (atds_mgwr) and regular multiscale GWR (multiscale_gwr)
* Introducing generalized GWR for binomial (bionomial and quasibinomial families).
* Introducing GAM/GWR with gradient descent boosting.
* Improving kernels for spatio-temporal GWR to introduce seasonnality (GDT)
* Removing unused experimental General kernel Product functions (GDX,GDC)
* Correcting a bug for MGWRSAR_x_x_x models with spatial autocorrelation when SE=TRUE.
* Correcting a bug for GWR model with a single explanatory variable when doMC=TRUE.
* Correcting a bug for formula without explicit names of variable (like "~.").
* Correcting a bug in GWR with parallel computation (split error in gwr_beta)


# 1.0.5 2023-11-16
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


# 1.0.4 2023-03-01
* Improved kernel_matW function by reintroducing possibility of island with no neighbours for non adaptie kernel
* Introduced an experimental version of a backfitting algorithm for the estimation of the Multiscale GWR with a selection of the bandiwth of each variable by cross validation (LOOCV). The predict_mgwrsar function allows to make predictions from a model of the 'multiscale_gwr' class.

# 1.0.3 2020-11-18
* Fixed a bug for plot_mgwrsar
* Fixed a bug for computing SE and edf when system is computationally singular
* Fixed a bug for the function plot_mgwrsar coming from the way of naming the columns of the coordinates matrix returned in the MGWRSAR function.
* Added a warning for cases where isgcv=TRUE and SE=TRUE when calling the MGWRSAR function.
* Improvement of plot_effects function to plot also effect of spatially non-varying parameters
* Improvement of plot_mgwrsar function to allow control of leaflet tile and number of quantile in legend.
* WARNING:  function to predict with target points using TP options have to be checked for model with spatial autocorrelation.


# 1.0.2 2020-06-01

* Fixed a bug in the kernel_matW function when there are duplicted spatal coordinates: duplicated coordinates are jittered with a warning if this is the case.

# 1.0.1 2020-05-04

* All models of mgwrsar package based on local linear regression can now be estimated using a target points set. Several functions that allows to choose an optimal set of target points to obtain a faster approximation of GWR coefficients has been added.

* Predictions on new data can now be done using the jacknife estimation method instead of spatial extrapolation of local coefficients from a preliminarily estimated model. Only the optimal value of the bandwidth modeled from the initial data is then used. In the function 'predict_mgwrsar', if the parameter method_pred = 'TP' (default), the prediction is done by recalculating a MGWRSAR model with the new data as target points keeping the bandwidth at the optimal value chosen with the training data, otherwise if method_pred= ('tWtp_model', 'model', 'shepard') then a matrix is used for the spatial extrapolation of the estimated coefficients, and prediction are done using these extrapolated coefficients (as in the previous version of mgwrsar package).

* 'KNN' function is deprecated and replaced by 'kernel_matW' function that allows to build spatial weight matrix and interaction matrix based on General Kernel Product. In kernel_matW function it's possible to specify the maximum number of neighbors to consider in gaussian kernel (rough gaussian kernel) to increase speed and sparsity of weights matrix.

* Fast computation of local OLS coefficients in the previous version of mgwrsar package (0.1) uses non pivotal computation that may provides undesirable results in presence of strong colinearity. The RCCP 'fastlmLLT_C' function has been replaced by R native lm.fit function in this realease.

* A new ploting function has been added: plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables

# 0.1 2018-05-11
* First release.

