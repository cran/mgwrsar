if(FALSE){

library(kableExtra)

i=1
mymodels=list(model_GWR_TP16f,model_GWR_TP8f,model_GWR_TP4f,model_GWR_fullf)
for(j in 1:length(mymodels)){
model<-mymodels[[j]]
n<-length(model$Y)
Beta_mean<-apply(model$Betav, 2,function(x) round(mean(x,na.rm=TRUE),2))
Beta_range <- apply(model$Betav, 2, function(x) paste0('[', paste(c(round(min(x),2),round(max(x),2)),collapse=' : '),']'))
Beta<-t(rbind(Beta_mean,Beta_range))
Beta<-Beta %>% melt() %>% arrange(X1) %>% select(value)
#names(Beta)<-'$\\beta$'

t_mean<-apply(model$Betav/model$sev, 2,function(x) round(mean(x,na.rm=TRUE),2))
t_range <- apply(model$Betav/model$sev, 2, function(x) paste0('[', paste(c(round(min(x,na.rm=TRUE),2),round(max(x,na.rm=TRUE),2)),collapse=' : '),']'))
tt<-t(rbind(t_mean,t_range))
tt<-tt %>% melt() %>% arrange(X1) %>% select(X1,value)
myrownamesV <-tt %>% select(X1)
myrownamesV<-as.character(myrownamesV$X1)
myrownamesV[seq(2,length(myrownamesV),by=2)]<-''
tv<-tt %>% select(value)
#names(tv)<-'Student t'
betav<-cbind(Beta,tv)
betac<-round(model$Betac,2)
tt<-round(model$Betac/model$se,2)
betac<-data.frame(betac,tt)
myrownames<-c(myrownamesV,row.names(betac))

names(betav)<-names(betac)<-c('$\\beta$','Student t')

globMod<-rbind(
  c('n :',length(model$Y)),
  c("Target Points:", ifelse(is.null(model$TP) | length(model$TP)==n,'NO',paste0('YES ',length(model$TP),' / ',n))),
  c("Edf: ", round(model$edf,2)),
  c('Kernel: ', model$kernels),
  c('adaptive: ', ifelse(model$adaptive,'YES','NO')),
  c('Rough kernel:',ifelse(model$NN<n & (!model$adaptive | (model$adaptive & model$kernels=='gauss')) ,paste0('YES, ',model$NN,' nb'),'NO')),
  c('Bandwidth: ', model$H),
  c('Computation time: ', paste0(round(model$ctime,2),'s')),
  c("AIC: ", round(2*n*log(sqrt(model$SSR/n)) + n*log(2*pi) + n + model$tS,2)),
  c("In sample RMSE: ",model$RMSEn),
  c("Out sample RMSPE: ",'XXXX')
)
globMod<-cbind(globMod,rep('',nrow(globMod)))

mytab<-cbind(myrownames,rbind(betav,betac))
colnames(globMod)<-colnames(mytab)
mytab<-rbind(mytab,globMod)

if(i==1) {
  mytabF<-mytab
  } else mytabF<-cbind(mytabF,mytab[,-1])
i=i+1
}


nmodel=(ncol(mytabF)-1)/2
myheader<-c(" ", "GWR TP16" = 2,"GWR TP8" = 2,"GWR TP4" = 2,"GWR full" = 2)

mytab1<-kbl(mytabF,row.names =FALSE,linesep = "",align=c('l',rep('c',2*nmodel)), col.names =c('',rep(c('$\\beta$','Student t'),nmodel)),'latex') %>% add_header_above(myheader)%>% kable_styling(latex_options = c("repeat_header")) %>% pack_rows("Varying Coefficient", 1, 14) %>% pack_rows("non Varying Coefficient", 15, 19) %>% kable_styling(full_width = FALSE) %>% row_spec(seq(2,14,by=2), italic=T) %>% row_spec(20:30, align = 'l')
}


