## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)



## ----load_data----------------------------------------------------------------
library(mgwrsar)
## loading data example
data(mydata)
coord=as.matrix(mydata[,c("x_lat","y_lon")])

## ----GWR_NN-------------------------------------------------------------------

## without rough gaussian kernel
ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE))
(proc.time()-ptm1)[3]

## with rough gaussian kernel
ptm1<-proc.time()
model_GWR_grk<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE,NN=300))
(proc.time()-ptm1)[3]

summary(model_GWR$Betav)
summary(model_GWR_grk$Betav)


## ----GWR_TP-------------------------------------------------------------------
TP=find_TP(formula = 'Y_gwr~X1+X2+X3', data =mydata,coord=coord,K=6,type='residuals')
# only 60 targets points are used
length(TP)

ptm1<-proc.time()
model_GWR_tp<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE,TP=TP,kWtp=12))
(proc.time()-ptm1)[3]

ptm1<-proc.time()
model_GWR_tp_NN<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE,TP=TP,kWtp=12,NN=300))
(proc.time()-ptm1)[3]

summary(model_GWR$Betav)
summary(model_GWR_tp$Betav)
summary(model_GWR_tp_NN$Betav)


## ----Prediction1--------------------------------------------------------------

length_out=800
index_in=sample(1:1000,length_out)
index_out=(1:1000)[-index_in]

coord_in=coord[index_in,]
data_in=mydata[index_in,]

TP=find_TP(formula = 'Y_gwr~X1+X2+X3', data =data_in,coord=coord_in,K=6,type='residuals')
# only 60 targets points are used
length(TP)

model_GWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars=NULL,kernels=c('gauss'),H=10, Model = 'GWR',control=list(adaptive=T,TP=TP))
summary_mgwrsar(model_GWR_insample)

newdata=mydata[index_out,]
newdata_coord=coord[index_out,]
newdata$Y_mgwrsar_1_0_kv=0

Y_pred=predict_mgwrsar(model_GWR_insample, newdata=newdata, newdata_coord=newdata_coord,method_pred='tWtp_model')
head(Y_pred)
head(mydata$Y_gwr[index_out])
sqrt(mean((mydata$Y_gwr[index_out]-Y_pred)^2)) # RMSE

## ----Prediction2--------------------------------------------------------------

length_out=200
index_in=sample(1:1000,length_out)
index_out=(1:1000)[-index_in]

### Global Spatial Weight matrix W should be ordered by in_ sample (S) then out_sample
W=kernel_matW(H=4,kernels='rectangle',coord_i=rbind(coord[index_in,],coord[index_out,]),NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)

W_in=kernel_matW(H=4,kernels='rectangle',coord_i=coord[index_in,],NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)

model_MGWRSAR_1_0_kv_insample<-MGWRSAR(formula = 'Y_mgwrsar_1_0_kv~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars=NULL,kernels=c('gauss'),H=11, Model = 'MGWRSAR_1_0_kv',control=list(W=W_in,adaptive=TRUE,isgcv=F))
model_MGWRSAR_1_0_kv_insample$RMSE
summary_mgwrsar(model_MGWRSAR_1_0_kv_insample)

## without Best Linear Unbiased Predictor
newdata=mydata[index_out,]
newdata_coord=coord[index_out,]
newdata$Y_mgwrsar_1_0_kv=0

Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='YTC')
head(Y_pred)
RMSE_YTC=sqrt(mean((mydata$Y_mgwrsar_1_0_kv[index_out]-Y_pred)^2))
RMSE_YTC

## Using Best Linear Unbiased Predictor
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='BPN')
head(Y_pred)
RMSE_BPN=sqrt(mean((mydata$Y_mgwrsar_1_0_kv[index_out]-Y_pred)^2))
RMSE_BPN


## ----bandwidths_mgwrsar with TP and NN, eval=FALSE----------------------------
#  
#  ptm1<-proc.time()
#  mytab_TP_NN<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,Models=c('GWR'),Kernels=c('gauss'),control=list(TP=TP,NN=300,adaptive=FALSE),control_search=list())
#  (proc.time()-ptm1)[3]
#  
#  names(mytab_TP_NN)
#  names(mytab_TP_NN[['GWR_gauss']])
#  mytab_TP_NN[['GWR_gauss']]$config_model
#  
#  
#  ptm1<-proc.time()
#  mytab<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,Models=c('GWR'),Kernels=c('gauss'),control=list(adaptive=FALSE),control_search=list())
#  (proc.time()-ptm1)[3]
#  
#  names(mytab)
#  names(mytab[['GWR_gauss']])
#  mytab[['GWR_gauss']]$config_model
#  

