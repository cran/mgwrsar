## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data----------------------------------------------------------------
library(mgwrsar)
library(dplyr)
## loading data example
data(mydata)
coord=as.matrix(mydata[,c("x_lat","y_lon")])
## Creating a spatial weight matrix (sparce dgCMatrix) of 8 nearest neighbors with 0 in diagonal
W=kernel_matW(H=4,kernels='rectangle',coord_i=coord,NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)

## ----GWR1---------------------------------------------------------------------

### without parallel computing with distance computation for all points
  ptm1<-proc.time()
  model_GWR0<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.13, Model = 'GWR',control=list(SE=T))
  (proc.time()-ptm1)[3]
  
  summary_mgwrsar(model_GWR0)


### with parallel computing with distance computation for all points
ptm1<-proc.time()
model_GWR0_parallel<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=T,ncore=2,doMC=TRUE))
(proc.time()-ptm1)[3]

# W<-spdep::mat2listw(W)
# 
# gp2mM <- lagsarlm(Y_gwr ~ X1+X2+X3,data=mydata, lW, method="Matrix")
# 
# mdo<-s2sls(Y_gwr ~ X1+X2+X3,data=mydata,Produc,W)
# 
# mm<-lagmess(Y_gwr ~ X1+X2+X3,data=mydata, lW)

## how to speed up computation using rough kernels (0 weights for nearest neighbors > 300th position)
ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE,NN=300))
(proc.time()-ptm1)[3]

summary_mgwrsar(model_GWR)

## how to speed up computation using rough kernels ( 0 weights for nearest neighbors > 300th position) 
TP=find_TP(formula = 'Y_gwr~X1+X2+X3', data =mydata,coord=coord,K=10,type='residuals')
# only 90 targets points are used
length(TP)

## how to speed up computation using Target points
ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(SE=TRUE,TP=TP))
(proc.time()-ptm1)[3]

## how to speed up computation using rough kernels ( 0 weights for nearest neighbors > 300th position)  and Target points

ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.1, Model = 'GWR',control=list(SE=TRUE,NN=300,TP=TP))
(proc.time()-ptm1)[3]

## how to speed up computation using adapative kernels 
ptm1<-proc.time()
model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('bisq'),H=20, Model = 'GWR',control=list(SE=TRUE,adaptive=TRUE))
(proc.time()-ptm1)[3]

## how to speed up computation using adapative kernels and Target points
  ptm1<-proc.time()
  model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('bisq'),H=20, Model = 'GWR',control=list(SE=TRUE,adaptive=TRUE,TP=TP))
  (proc.time()-ptm1)[3]


library(microbenchmark)
res=microbenchmark(MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.1, Model = 'GWR',control=list(SE=TRUE)),MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.1, Model = 'GWR',control=list(SE=TRUE,NN=300,TP=TP)),times=2)
res

res=microbenchmark(MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.1, Model = 'GWR',control=list(SE=TRUE)),MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('bisq'),H=20, Model = 'GWR',control=list(SE=TRUE,adaptive=T,TP=TP)),times=2)
res

## ----plot1--------------------------------------------------------------------
summary_mgwrsar(model_GWR0)
plot_mgwrsar(model_GWR0,type='B_coef',var='X2')
plot_mgwrsar(model_GWR0,type='t_coef',var='X2')

## ----plot2--------------------------------------------------------------------
plot_effect(model_GWR0,title='Effects')

## ----spgwr, eval=FALSE--------------------------------------------------------
#  library(spgwr)
#  mydataSP=mydata
#  coordinates(mydataSP)=c('x_lat','y_lon')
#  ptm1<-proc.time()
#  model_spgwr <- gwr(Y_gwr~X1+X2+X3, data=mydataSP, bandwidth=0.13,hatmatrix=TRUE)
#  (proc.time()-ptm1)[3]
#  
#  head(model_spgwr$SDF@data$gwr.e)
#  model_spgwr
#  
#  all(abs(model_GWR0$residuals-model_spgwr$SDF@data$gwr.e)<0.00000000001)
#  

## ----GWR2---------------------------------------------------------------------
model_GWR_loo_no_outlier<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('bisq'),H=20, Model = 'GWR',control=list(isgcv=TRUE,adaptive=TRUE,remove_local_outlier=TRUE,outv=0.01))
summary_mgwrsar(model_GWR_loo_no_outlier)

### leave-one out CV estimation
model_GWR_loo<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=20, Model = 'GWR',control=list(isgcv=TRUE,adaptive=T))
summary_mgwrsar(model_GWR_loo)


## ----MGWR1--------------------------------------------------------------------

model_MGWR<-MGWRSAR(formula = 'Y_mgwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars='X3',kernels=c('gauss'),H=20, Model = 'MGWR',control=list(SE=T,adaptive=T))
summary_mgwrsar(model_MGWR)
             

## ----mgwrsar_0_kc_kv----------------------------------------------------------
mgwrsar_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars='X2',kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_kc_kv',control=list(SE=F,adaptive=T,W=W))
summary_mgwrsar(mgwrsar_0_kc_kv)

## ----mgwrsar_0_0_kv-----------------------------------------------------------
mgwrsar_0_0_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_0_kv',control=list(SE=F,adaptive=T,W=W))
summary_mgwrsar(mgwrsar_0_0_kv)

## ----mgwrsar_1_0_kv-----------------------------------------------------------
mgwrsar_1_0_kv<-MGWRSAR(formula = 'Y_mgwrsar_1_0_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=20, Model = 'MGWRSAR_1_0_kv',control=list(SE=F,adaptive=T,W=W))
summary_mgwrsar(mgwrsar_1_0_kv)

## ----mgwrsar_1_kc_kv----------------------------------------------------------
mgwrsar_1_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_1_kc_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars='X2',kernels=c('gauss'),H=20, Model = 'MGWRSAR_1_kc_kv',control=list(SE=F,adaptive=T,W=W))
summary_mgwrsar(mgwrsar_1_kc_kv)

## ----mgwrsar_1_kc_0-----------------------------------------------------------
mgwrsar_1_kc_0<-MGWRSAR(formula = 'Y_mgwrsar_1_kc_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars=c('Intercept','X1','X2','X3'),kernels=c('gauss'),H=20, Model = 'MGWRSAR_1_kc_0',control=list(SE=F,adaptive=T,W=W))
summary_mgwrsar(mgwrsar_1_kc_0)

## ----Prediction1--------------------------------------------------------------
length_out=800
index_in=sample(1:1000,length_out)
index_out=(1:1000)[-index_in]

model_GWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars=NULL,kernels=c('gauss'),H=8, Model = 'GWR',control=list(adaptive=T))
summary_mgwrsar(model_GWR_insample)

newdata=mydata[index_out,]
newdata_coord=coord[index_out,]
newdata$Y_mgwrsar_1_0_kv=0

Y_pred=predict_mgwrsar(model_GWR_insample, newdata=newdata, newdata_coord=newdata_coord)
head(Y_pred)
head(mydata$Y_gwr[index_out])
sqrt(mean((mydata$Y_gwr[index_out]-Y_pred)^2)) # RMSE

model_MGWR_insample<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars='X3',kernels=c('gauss'),H=8, Model = 'MGWR',control=list(adaptive=T))
summary_mgwrsar(model_MGWR_insample)

newdata=mydata[index_out,]
newdata_coord=coord[index_out,]
newdata$Y_mgwrsar_1_0_kv=0

## Prediction with method_pred='tWtp' 
Y_pred=predict_mgwrsar(model_MGWR_insample, newdata=newdata, newdata_coord=newdata_coord)
head(Y_pred)
head(mydata$Y_gwr[index_out])
sqrt(mean((mydata$Y_gwr[index_out]-Y_pred)^2)) # RMSE


## ----Prediction2--------------------------------------------------------------

length_out=800
index_in=sample(1:1000,length_out)
index_out=(1:1000)[-index_in]

### Global Spatial Weight matrix W should be ordered by in_sample (S) then out_sample
W_in_out=kernel_matW(H=4,kernels='rectangle',coord_i=rbind(coord[index_in,],coord[index_out,]),NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)
W_in=W_in_out[1:length(index_in),1:length(index_in)]
W_in=mgwrsar::normW(W_in)
newdata=mydata[index_out,]
newdata_coord=coord[index_out,]


### SAR Model
model_SAR_insample<-MGWRSAR(formula = 'Y_sar~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars=NULL,kernels=NULL,H=11, Model = 'SAR',control=list(W=W_in))
model_SAR_insample$RMSE
summary_mgwrsar(model_SAR_insample)

## without Best Linear Unbiased Predictor
# prediction YTC
Y_pred=predict_mgwrsar(model_SAR_insample, newdata=newdata, newdata_coord=newdata_coord,W=W_in_out,type='YTC')
head(Y_pred)
RMSE_YTC=sqrt(mean((mydata$Y_sar[index_out]-Y_pred)^2))
RMSE_YTC

## Using Best Linear Unbiased Predictor
Y_pred=predict_mgwrsar(model_SAR_insample, newdata=newdata, newdata_coord=newdata_coord,W=W_in_out,type='BPN')
head(Y_pred)
RMSE_BPN=sqrt(mean((mydata$Y_sar[index_out]-Y_pred)^2))
RMSE_BPN


#### MGWRSAR_1_0_kv
model_MGWRSAR_1_0_kv_insample<-MGWRSAR(formula = 'Y_mgwrsar_1_0_kv~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars=NULL,kernels=c('gauss'),H=11, Model = 'MGWRSAR_1_0_kv',control=list(W=W_in,adaptive=TRUE,isgcv=F))
model_MGWRSAR_1_0_kv_insample$RMSE
summary_mgwrsar(model_MGWRSAR_1_0_kv_insample)



# prediction with method_pred='tWtp'
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W_in_out,type='BPN',method_pred='tWtp')
head(Y_pred)
RMSE_tWtp_BPN=sqrt(mean((mydata$Y_mgwrsar_1_0_kv[index_out]-Y_pred)^2))
RMSE_tWtp_BPN


## Using Best Linear Unbiased Predictor with method_pred='model' 
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='BPN',method_pred='model' )
head(Y_pred)
RMSE_model_BPN=sqrt(mean((mydata$Y_mgwrsar_1_0_kv[index_out]-Y_pred)^2))
RMSE_model_BPN

## Using Best Linear Unbiased Predictor with method_pred='sheppard' 
W_out=kernel_matW(H=4,kernels='rectangle',coord_i=coord[index_out,],NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)
Y_pred=predict_mgwrsar(model_MGWRSAR_1_0_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='BPN',method_pred='sheppard',k_extra=8)
head(Y_pred)
RMSE_sheppard_BPN=sqrt(mean((mydata$Y_mgwrsar_1_0_kv[index_out]-Y_pred)^2))
RMSE_sheppard_BPN



## ----Prediction3--------------------------------------------------------------

length_out=800
index_in=sample(1:1000,length_out)
index_out=(1:1000)[-index_in]

### Global Spatial Weight matrix W should be ordered by in_ sample (S) then out_sample
W=kernel_matW(H=4,kernels='rectangle',coord_i=rbind(coord[index_in,],coord[index_out,]),NN=4,adaptive=TRUE,diagnull=TRUE,rowNorm=T)

W_in=W[index_in,index_in]
W_in=mgwrsar::normW(W_in)

model_MGWRSAR_1_kc_kv_insample<-MGWRSAR(formula = 'Y_mgwrsar_1_kc_kv~X1+X2+X3', data = mydata[index_in,],coord=coord[index_in,], fixed_vars='X3',kernels=c('gauss'),H=11, Model = 'MGWRSAR_1_kc_kv',control=list(W=W_in,adaptive=TRUE,isgcv=F))
model_MGWRSAR_1_kc_kv_insample$RMSE
summary_mgwrsar(model_MGWRSAR_1_kc_kv_insample)

## without Best Linear Unbiased Predictor
newdata=mydata[index_out,]
newdata_coord=coord[index_out,]
newdata$Y_mgwrsar_1_kc_kv=0

Y_pred=predict_mgwrsar(model_MGWRSAR_1_kc_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='YTC')
head(Y_pred)
RMSE_YTC=sqrt(mean((mydata$Y_mgwrsar_1_kc_kv[index_out]-Y_pred)^2))
RMSE_YTC

## Using Best Linear Unbiased Predictor
Y_pred=predict_mgwrsar(model_MGWRSAR_1_kc_kv_insample, newdata=newdata, newdata_coord=newdata_coord,W=W,type='BPN')#,method_pred='sheppard')
head(Y_pred)
RMSE_BPN=sqrt(mean((mydata$Y_mgwrsar_1_kc_kv[index_out]-Y_pred)^2))
RMSE_BPN


## ----bandwidths_search,eval=FALSE---------------------------------------------
#  
#  ######################
#  #### Finding bandwidth by hand
#  #####################
#  model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=0.03, Model = 'GWR',control=list(isgcv=TRUE,adaptive=FALSE))
#  summary_mgwrsar(model_GWR)
#  
#    myCV<-function(H){
#      model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H, Model = 'GWR',control=list(isgcv=TRUE,adaptive=FALSE))
#      cat('... Try h=',H,' ')
#      CV=model_GWR$RMSE
#      CV
#    }
#  
#  resCV=optimize(myCV,lower=0.01,upper=0.4)
#  
#  resCV
#  
#  ## model with optimal bandwith with gaussian kernel
#  
#  model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss'),H=resCV$minimum, Model = 'GWR',control=list(isgcv=F,adaptive=F))
#  summary_mgwrsar(model_GWR)

## ----bandwidths_mgwrsar1,eval=FALSE-------------------------------------------
#  mytab<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=c('Intercept','X1'),Models=c('GWR','MGWR'),candidates_Kernels=c('bisq','gauss'),control=list(NN=300,adaptive=TRUE),control_search=list())
#  
#  names(mytab)
#  names(mytab[['GWR_bisq_adaptive']])
#  
#  mytab[['GWR_bisq_adaptive']]$config_model
#  mytab[['GWR_bisq_adaptive']]$CV
#  summary(mytab[['GWR_bisq_adaptive']]$model$Betav)
#  
#  mybestmodel=mytab[['GWR_gauss_adaptive']]$model

## ----bandwidths_mgwrsar2, eval=FALSE------------------------------------------
#  mytab2<-bandwidths_mgwrsar(formula='Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,Models=c('MGWRSAR_0_0_kv'),candidates_Kernels=c('bisq'),control=list(adaptive=TRUE,NN=500),control_search=list(search_W=TRUE,kernel_w='rectangle',adaptive_W=TRUE))
#  
#  names(mytab2)
#  mytab2[['MGWRSAR_0_0_kv_bisq_adaptive']]$config_model
#  mytab2[['MGWRSAR_0_0_kv_bisq_adaptive']]$CV
#  mytab2[['MGWRSAR_0_0_kv_bisq_adaptive']]$model$Betac
#  summary(mytab2[['MGWRSAR_0_0_kv_bisq_adaptive']]$model$Betav)

## ----EXPERIMENTAL, eval=TRUE--------------------------------------------------
## space + time kernel
mytime_index=sort(rep(1:500,2))
mytime_index[1:150]
W=kernel_matW(H=8,kernels='rectangle',coord_i=coord,NN=10,adaptive=TRUE,diagnull=TRUE,rowNorm=T)

  
model_MGWRSART_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,coord=coord, fixed_vars=c('Intercept'),kernels=c('gauss','gauss'),H=c(50,50), Model = 'MGWRSAR_0_kc_kv',control=list(Z=mytime_index,W=W,adaptive=c(TRUE,TRUE),Type='GDT'))

summary_mgwrsar(model_MGWRSART_0_kc_kv)


### space  + continuous variable kernel
model_GWRX<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord, fixed_vars=NULL,kernels=c('gauss','gauss'),H=c(120,120), Model = 'GWR',control=list(Z=mydata$X2,Type='GDX',adaptive=c(T,T)))
summary_mgwrsar(model_GWRX)


### space + categorical kernel (Li and Racine 2010)
Z=1+as.numeric(mydata$X2>quantile(mydata$X2,0.9))*2+as.numeric(mydata$X2<=quantile(mydata$X2,0.1))
table(Z)

model_MGWRSARC_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X3', data = mydata,coord=coord, fixed_vars=c('Intercept'),kernels=c('gauss','li','li','li'),H=c(120,0.1,0.9,0.9), Model = 'MGWRSAR_0_kc_kv',control=list(Z=Z,W=W,Type='GDC',adaptive=c(T,F,F,F)))
summary_mgwrsar(model_MGWRSARC_0_kc_kv)

