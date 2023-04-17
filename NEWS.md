NEWS/ChangeLog
-----------------------------

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

* Predictions on new data can now be done using the jacknife estimation method instead of spatial extrapolation of local coefficients from a preliminarily estimated model. Only the optimal value of the bandwidth modeled from the initial data is then used. In the function 'predict_mgwrsar', if the parameter method_pred = 'TP' (default), the prediction is done by recalculating a MGWRSAR model with the new data as target points keeping the bandwidth at the optimal value chosen with the training data, otherwise if method_pred= ('tWtp_model', 'model', 'sheppard') then a matrix is used for the spatial extrapolation of the estimated coefficients, and prediction are done using these extrapolated coefficients (as in the previous version of mgwrsar package).

* 'KNN' function is deprecated and replaced by 'kernel_matW' function that allows to build spatial weight matrix and interaction matrix based on General Kernel Product. In kernel_matW function it's possible to specify the maximum number of neighbors to consider in gaussian kernel (rough gaussian kernel) to increase speed and sparsity of weights matrix.

* Fast computation of local OLS coefficients in the previous version of mgwrsar package (0.1) uses non pivotal computation that may provides undesirable results in presence of strong colinearity. The RCCP 'fastlmLLT_C' function has been replaced by R native lm.fit function in this realease.

* A new ploting function has been added: plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables

# 0.1 2018-05-11
* First release.

