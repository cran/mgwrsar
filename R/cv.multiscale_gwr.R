#' multiscale_gwr.cv
#' to be documented (experimental)
#' @usage multiscale_gwr.cv(dataName, argDataName="data", target='Y', K=5, regFun, par_model,
#' par_model2=NULL,regFun2=NULL, predFun, args_predNames, extra_args_pred=NULL,
#' namesXtraArgs2Split=NULL,myseed=1)
#' @param dataName  character, name of the data
#' @param argDataName character, generic name to use as data name.
#' @param target  character, name of variable to explain
#' @param K integer, number of folds for cross validation
#' @param regFun character, name of the estimation function
#' @param par_model named list with the arguments for the estimation function
#' @param predFun  character, name of the prediction function
#' @param args_predNames named list with the arguments for the prediction function
#' @param extra_args_pred named list with extra arguments for non generic prediction function
#' @param par_model2 to be documented
#' @param regFun2  to be documented
#' @param myseed seed for random number.
#' @param namesXtraArgs2Split character, names of the objects in extra_args_pred that need to be split for cross validation.
multiscale_gwr.cv=function(dataName, argDataName="data", target='Y', K=5, regFun, par_model,par_model2=NULL,regFun2=NULL, predFun, args_predNames, extra_args_pred=NULL, namesXtraArgs2Split=NULL,myseed=1){
  data=get(dataName)
  extra_args_pred_copy=extra_args_pred
  n=nrow(data)
  set.seed(myseed)
  cvFoldsList=createFolds(1:n, k=K, list=TRUE, returnTrain=FALSE)
  prediction=numeric(nrow(data))
  Resultat=list()

  for(i in 1:K){
    testset=data[cvFoldsList[[i]],]
    assign(dataName,data[-cvFoldsList[[i]],])
    param_model=eval(parse(text=par_model))

    if(!is.null(param_model$coord)){
      if(nrow(param_model$coord)==nrow(data)) param_model$coord=param_model$coord[-cvFoldsList[[i]],]
    }

    # Estimation

    #param_model[[argDataName]]=dataName
    # Model=do.call(regFun, args_model)
    Model=do.call(regFun, param_model)
    error=residuals(Model)
    if(regFun=='gwr.multiscale') {
      if(!is.null(par_model2)) param_model2=eval(parse(text=par_model2))
      error=Model$SDF$residual
      param_model2$coord<-param_model2$coord[-cvFoldsList[[i]],]
      param_model2$data<-param_model2$data[-cvFoldsList[[i]],]
      data2=data@data
      Model2=do.call(regFun2, param_model2)
      m=ncol(Model2$Betav)
      nx=colnames(Model2$Betav)
      Model2$Betav=as.matrix(Model$SDF@data[,1:m])
      colnames(Model2$Betav)=nx
      Model2$TP=1:nrow(Model2$Betav)
      Model=Model2
      Model$Model ='multiscale_GWR'
    } else data2=data


    # Prediction
    args_pred=list(model=Model,newdata=testset)
    names(args_pred)=args_predNames # Attention l'ordre est important
    if(!is.null(extra_args_pred)) {
      if(!is.null(namesXtraArgs2Split)){
        for(arg in namesXtraArgs2Split){
          if(ncol(extra_args_pred[[arg]]) > 1){
            extra_args_pred[[arg]]=extra_args_pred_copy[[arg]][cvFoldsList[[i]],]
          } else extra_args_pred[[arg]]=extra_args_pred_copy[[arg]][cvFoldsList[[i]]]
        }
      }
      args_pred=c(args_pred, extra_args_pred)
    }

    pred=do.call(predFun,args_pred)
    prediction[cvFoldsList[[i]]]=pred
    # sub_prediction=predict_mgwrsar(mygwr2,newdata = sub_testset, newdata_coord = coords_subTest)
    Resultat$Fold_index[[i]]=cvFoldsList[[i]]
    Resultat$Fold_Prediction[[i]]=pred
    # Resultat$RMSE_in[[i]]=Model$RMSE
    Resultat$RMSE_in[[i]]=sqrt(mean(error^2))
    Resultat$MAPE_in[[i]]=mean(abs(error/data2[-cvFoldsList[[i]],target]))*100
  }
  # pmse_model=sqrt(mean((testset[[target]]-prediction)^2))
  Resultat$Global$RMSE_insample=mean(unlist(Resultat$RMSE_in))
  Resultat$Global$MAPE_insample=mean(unlist(Resultat$MAPE_in))
  Resultat$Global$RMSE_outsample=sqrt(mean((data2[[target]]-prediction)^2))
  Resultat$Global$MAPE_outsample=mean(abs((data2[[target]]-prediction)/data[[target]]))*100
  return(Resultat)
}
