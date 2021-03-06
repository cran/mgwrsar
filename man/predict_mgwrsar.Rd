% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.mgwrsar.R
\name{predict_mgwrsar}
\alias{predict_mgwrsar}
\title{mgwrsar Model Predictions}
\usage{
predict_mgwrsar(model, newdata, newdata_coord, W = NULL,
type = "BPN", h_w = 100,kernel_w = "knn", k_extra = 12,
kernel_extra = "sheppard")
}
\arguments{
\item{model}{a model of mgwrsar class.}

\item{newdata}{a matrix or data.frame of new data.}

\item{newdata_coord}{a matrix of new coordinates.}

\item{W}{the spatial weight matrix for models with  spatial autocorrelation.}

\item{type}{Type for BLUP estimator, default "BPN". If NULL use predictions without spatial bias correction.}

\item{h_w}{bandwidth for constructing W, if W is NULL.}

\item{kernel_w}{kernel for constructing W, if W is NULL.}

\item{k_extra}{number of neighboors for local parameter extrapolation, default 12.}

\item{kernel_extra}{kernel for local parameter extrapolation, default sheppard kernel.}
}
\value{
A vector of predictions.
}
\description{
predict_mgwrsar is a function for computing predictions of a mgwrsar models. It uses Best Linear Unbiased Predictor for mgwrsar models with spatial autocorrelation.
}
\examples{
\donttest{
library(mgwrsar)
data(mydata)
coord=as.matrix(mydata[,c("x_lat","y_lon")])
W=KNN(coord,2)
model_GWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[1:800,],
coord=coord[1:800,],fixed_vars=NULL,kernels=c('gauss_adapt'),H=50,
 Model = 'GWR',control=list())
Y_pred=predict_mgwrsar(model_GWR_insample, newdata=mydata[801:1000,],
newdata_coord=coord[801:1000,],k_extra = 8, kernel_extra = "sheppard")
head(Y_pred)
head(mydata$Y_gwr[801:1000])
sqrt(mean((mydata$Y_gwr[801:1000]-Y_pred)^2))

## predict with spatial autocorrelation
model_MGWRSAR_1_0_kv_insample<-MGWRSAR(formula = 'Y_mgwrsar_1_0_kv~X1+X2+X3',data = mydata[1:800,],
coord=coord[1:800 ,], fixed_vars=NULL,kernels=c('gauss_adapt'),H=50,
Model = 'MGWRSAR_1_0_kv',control=list(W=W[1:800,1:800],Lambdacor=TRUE,SE=TRUE))
summary_mgwrsar(model_MGWRSAR_1_0_kv_insample)

## with BLUP
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=mydata[801:1000,],
newdata_coord=coord[801:1000,], k_extra = 12,  W = W,
 type = "BPN",kernel_extra = "sheppard")
head(Y_pred)
head(mydata$Y_gwr[801:1000])
sqrt(mean((mydata$Y_gwr[801:1000]-Y_pred)^2))

## without BLUP
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=mydata[801:1000,],
newdata_coord=coord[801:1000,], k_extra = 12,  W = W,
type = "TC",kernel_extra = "sheppard")
head(Y_pred)
head(mydata$Y_mgwrsar_1_0_kv[801:1000])
sqrt(mean((mydata$Y_mgwrsar_1_0_kv[801:1000]-Y_pred)^2))
}
}
\seealso{
MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, kernelW_C
}
