clean_colnames<-function(x){
  colnames(x)<-gsub('\\(','',colnames(x))
  colnames(x)<-gsub('\\^','',colnames(x))
  gsub('\\)','',colnames(x))
}
