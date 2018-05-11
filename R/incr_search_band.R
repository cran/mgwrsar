#' incr_search_band
#' to be documented
#' @usage incr_search_band(f,lower,discrete=FALSE,n,...)
#' @param f  to be documented
#' @param lower  to be documented
#' @param discrete  to be documented
#' @param n  to be documented
#' @param ...  to be documented
#' @keywords internal
#' @return to be documented
incr_search_band <-
function(f,lower,discrete=FALSE,n,...){
    fi<-function(x,pas=0.3,discrete=FALSE) ifelse(discrete,x+ceiling(x*pas),x+x*pas)
    f0<-f(lower,...)
    mymin=f0
    opt=TRUE
    onemin=FALSE
    lo=lower
    starting=lo
    i<-1
while(f(fi(lower,pas=0.3,discrete=discrete),...)<=mymin & ifelse(discrete,lower<trunc(n*0.65),TRUE)) { cat(lower,' . ');
mymin=f(fi(lower,pas=0.3,discrete=discrete),...);lo=lower;lower=fi(lower,pas=0.3,discrete=discrete);i=i+1;
}
    i<-1
    lower=lo
    mymin=f(lo,...)
    while(f(fi(lower,pas=0.1,discrete=discrete),...)<=mymin &ifelse(discrete,lower<trunc(n*0.65),TRUE)) {
    cat(lower,' . ');
  mymin=f(fi(lower,pas=0.1,discrete=discrete),...);lo=lower;lower=fi(lower,pas=0.1,discrete=discrete);i=i+1;
    }
     i<-1
    lower=lo
    mymin=f(lo,...)
    while(f(fi(lower,pas=0.01,discrete=discrete),...)<=mymin & ifelse(discrete,lower<trunc(n*0.65),TRUE)) {cat(lower,' . ');
    mymin=f(fi(lower,pas=0.01,discrete=discrete),...);lo=lower;lower=fi(lower,pas=0.01,discrete=discrete);i=i+1;
    #cat('i ',i,'  v',mymin,'  lower',lower,'  lo',lo,'\n')
    }
    if(!discrete){
        i<-1
        lower=lo
        mymin=f(lo,...)
        while(f(fi(lower,pas=0.02,discrete=discrete),...)<=mymin & ifelse(discrete,lower<trunc(n*0.65),TRUE)) { cat(lower,' . ');
        mymin=f(fi(lower,pas=0.02,discrete=discrete),...);lo=lower;lower=fi(lower,pas=0.02,discrete=discrete);i=i+1;
        #cat('i ',i,'  v',mymin,'  lower',lower,'  lo',lo,'\n')
        cat('.')
        }
    }
    list(minimum=lower,objective=mymin)
}
