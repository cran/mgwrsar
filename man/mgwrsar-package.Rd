\name{mgwrsar-package}
\alias{mgwrsar-package}
\alias{mgwrsar}
\docType{package}
\title{
\packageTitle{mgwrsar}
}
\description{
mgwrsar package proposes functions for estimating linear and local linear model with spatial autocorrelation. It allows to estimate linear and Spatial Autoregressive models with spatially varying coefficients. Models that mixed spatially varying and stationary coefficients can also be estimated.
\packageDESCRIPTION{mgwrsar}
\packageIndices{mgwrsar}
}
\details{
The DESCRIPTION file:
\packageDescription{mgwrsar}
}
\author{
Ghislain Geniaux \email{ghislain.geniaux hat inra.fr}
Davide Martinetti \email{davide.martinetti hat inra.fr}
}
\references{
Geniaux, G. and Martinetti, D.(2017). A new method for dealing simultaneously with spatial autocorrelation and spatial heterogeneity in regression models. Regional Science and Urban Economics. (https://doi.org/10.1016/j.regsciurbeco.2017.04.001)
}
\keyword{GWR, mixed GWR, LWR, local regression, SAR, spatial autocorrelation}
\seealso{
\code{\link[McSpatial]{lwr}}
\code{\link[locfit]{locfit}}
}
\examples{
library(mgwrsar)
## loading data example
data(data_mgwrsar)
coord=as.matrix(mydata[,c("x_lat","y_lon")])
## Creating a spatial weigth matrix (sparce dgCMatrix) of 8 nearest neighbors
W=KNN(coord,8)

ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model = 'GWR',
control=list(SE=TRUE,doMC=FALSE))
(proc.time()-ptm1)[3]

summary_mgwrsar(model_GWR)
plot_mgwrsar(model_GWR,type='B_coef',var='X2')
plot_mgwrsar(model_GWR,type='t_coef',var='X2')
}
