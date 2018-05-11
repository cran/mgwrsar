
// includes from the plugin
#include <RcppEigen.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;
typedef Map<VectorXd> MapVecd;
inline MatrixXd AtA(const MatrixXd& A) {
int n(A.cols());
return MatrixXd(n,n).setZero().selfadjointView<Lower>()
.rankUpdate(A.adjoint());
}

// declarations
extern "C" {
SEXP Proj_C( SEXP HH, SEXP XX) ;
SEXP D_dense_C( SEXP xxi, SEXP yyi, SEXP xx, SEXP yy) ;
SEXP Sl_C( SEXP llambda, SEXP WW, SEXP iinv, SEXP aapprox) ;
SEXP INST_C( SEXP XX, SEXP WW, SEXP withlambda, SEXP llambda) ;
SEXP bisq_C( SEXP dd, SEXP hh, SEXP Minv) ;
SEXP bisq_knn_C( SEXP dd, SEXP hh) ;
SEXP bin_C( SEXP dd, SEXP hh, SEXP Minv) ;
SEXP gauss_C( SEXP dd, SEXP hh, SEXP Minv) ;
SEXP gauss_knn_C( SEXP dd, SEXP hh) ;
SEXP gauss_adapt_C( SEXP dd, SEXP hh) ;
SEXP Dx_dense_C( SEXP xxi, SEXP xx, SEXP TIME) ;
SEXP set0_conditionXD( SEXP XX, SEXP ZZ) ;
SEXP kernelW_C( SEXP XX, SEXP hh, SEXP MykernelS, SEXP isgcv_, SEXP Type, SEXP Minv, SEXP maxknn_, SEXP NmaxDist_, SEXP TIME, SEXP Decay, SEXP DDiagNull) ;
SEXP kernel_C( SEXP XX, SEXP J, SEXP hh, SEXP Mykernel, SEXP Minv, SEXP TIME, SEXP Decay, SEXP DDiagNull, SEXP normWW) ;
SEXP PhWY_C( SEXP YY, SEXP XX, SEXP WW, SEXP Wi) ;
SEXP fastlmLLT_C( SEXP XX, SEXP YY, SEXP SE_) ;
SEXP QRcpp_C( SEXP AA, SEXP bb) ;
SEXP QRcpp2_C( SEXP AA, SEXP bb, SEXP cc) ;
SEXP ApproxiW( SEXP WW, SEXP la, SEXP order) ;
SEXP mod( SEXP YY, SEXP XX, SEXP WW, SEXP XZZ, SEXP YZZ, SEXP Wi, SEXP LocalInst, SEXP ismethodB2SLS, SEXP ismethodMGWRSAR_1_kc_0, SEXP SE_) ;
}

static const R_CallMethodDef CallEntries[] = {
  {"ApproxiW",         (DL_FUNC) &ApproxiW,          3},
  {"bin_C",            (DL_FUNC) &bin_C,             3},
  {"bisq_C",           (DL_FUNC) &bisq_C,            3},
  {"bisq_knn_C",       (DL_FUNC) &bisq_knn_C,        2},
  {"D_dense_C",        (DL_FUNC) &D_dense_C,         4},
  {"Dx_dense_C",       (DL_FUNC) &Dx_dense_C,        3},
  {"fastlmLLT_C",      (DL_FUNC) &fastlmLLT_C,       3},
  {"gauss_adapt_C",    (DL_FUNC) &gauss_adapt_C,     2},
  {"gauss_C",          (DL_FUNC) &gauss_C,           3},
  {"gauss_knn_C",      (DL_FUNC) &gauss_knn_C,       2},
  {"INST_C",           (DL_FUNC) &INST_C,            4},
  {"kernel_C",         (DL_FUNC) &kernel_C,          9},
  {"kernelW_C",        (DL_FUNC) &kernelW_C,        11},
  {"mod",              (DL_FUNC) &mod,              10},
  {"PhWY_C",           (DL_FUNC) &PhWY_C,            4},
  {"Proj_C",           (DL_FUNC) &Proj_C,            2},
  {"QRcpp_C",          (DL_FUNC) &QRcpp_C,           2},
  {"QRcpp2_C",         (DL_FUNC) &QRcpp2_C,          3},
  {"set0_conditionXD", (DL_FUNC) &set0_conditionXD,  2},
  {"Sl_C",             (DL_FUNC) &Sl_C,              4},
  {NULL, NULL, 0}
};

void R_init_mgwrsar(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
// definition

SEXP Proj_C( SEXP HH, SEXP XX ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
MatrixXd X=Rcpp::as<MatrixXd>(XX);
MatrixXd H=Rcpp::as<MatrixXd>(HH);
MatrixXd res = H  * (H.adjoint()*H).inverse() * (H.transpose()* X);
return wrap(res);
END_RCPP
}


SEXP D_dense_C( SEXP xxi, SEXP yyi, SEXP xx, SEXP yy ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
VectorXd x=Rcpp::as<VectorXd>(xx);
VectorXd y=Rcpp::as<VectorXd>(yy);
const double xi= Rcpp::as<double>(xxi) ;
const double yi= Rcpp::as<double>(yyi) ;
VectorXd x1= x.array()-xi;
VectorXd y1= y.array()-yi;
VectorXd res= x1.array().square()+y1.array().square();
res = res.array().sqrt();
return wrap(res);
END_RCPP
}


SEXP Sl_C( SEXP llambda, SEXP WW, SEXP iinv, SEXP aapprox ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::SparseMatrix;

SparseMatrix<double>  W=Rcpp::as<SparseMatrix<double> >(WW);
double lambda =  as<double>(llambda);
bool inv = as<bool>(iinv);
bool approx = as<bool>(aapprox);
int n = W.rows();

SparseMatrix<double> I(n,n);
I.setIdentity();

Environment matr("package:Matrix");
Function solve = matr["solve"];

SparseMatrix<double> SW = I - lambda*W;

if(inv==true) {
    if(approx==true) {
                      SparseMatrix<double> W2 = W*W;
                      double lambda2 = lambda*lambda;
                      SparseMatrix<double> res = I + lambda*W + lambda2*W2 + lambda2*lambda*W*W2 + lambda2*lambda2*W2*W2;
                      return wrap(res);
                      } else {
                      SEXP res = solve(SW);
                      return wrap(res);
                      }
              } else {
                    SparseMatrix<double> res = SW;
                    return wrap(res);
}
END_RCPP
}


SEXP INST_C( SEXP XX, SEXP WW, SEXP withlambda, SEXP llambda ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

MatrixXd x=Rcpp::as<MatrixXd>(XX);
SparseMatrix<double>  W=Rcpp::as<SparseMatrix<double> >(WW);

double lambda =Rcpp::as<double>(llambda);
bool wl = Rcpp::as<bool>(withlambda);

Function ReorderX("int_prems");

SEXP xx=ReorderX(x);
MatrixXd X=as<MatrixXd>(xx);
int n = W.rows();
int m = X.cols();

Function Sl_C("Sl_C");
Function sd("sd");

Function cbind("cbind");
Environment Mat("package:Matrix");
Function cBind = Mat["cBind"];

double c0sum = X.col(0).sum();
SEXP sd0= sd(X.col(0));
double sdc0 = as<double>(sd0);
bool check = c0sum==n && sdc0==0;

SEXP iWW;
SEXP H;
MatrixXd iWX;
MatrixXd WX;
MatrixXd WWX;
MatrixXd WWWX;
MatrixXd HH;
if (wl) {iWW = Sl_C(lambda,W,true,true);
        SparseMatrix<double > iW = as<SparseMatrix<double> >(iWW);
        if (check)  {
        iWX = W*iW * (X.rightCols(m-1));
        } else {iWX = W* iW * X; }
        H = cbind(X,iWX);} else {
        WX = W*(X.rightCols(m-1));
        WWX = W*WX;
        WWWX = W*WWX;
        H = cBind(X,WX,WWX,WWWX);
        }
HH=as<MatrixXd>(H);
return wrap(HH);
END_RCPP
}


SEXP bisq_C( SEXP dd, SEXP hh, SEXP Minv ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;
VectorXd d=Rcpp::as<VectorXd>(dd);
double h = Rcpp::as<double>(hh);
int minv = Rcpp::as<int>(Minv);
Function sort("sort");
SEXP rkk =sort(d);
VectorXd zz=as<VectorXd>(rkk);
if(h<zz(minv)) h=zz(minv);

double c =15./16;
VectorXd x = d/h;
VectorXd y = x.cwiseAbs();
int n = y.size();
VectorXd z = c*(1-x.array().square()).array().square();
for(int j = 0; j < n; ++j)
{
if(y(j)>=1) z(j) = 0;
}
return wrap(z);
END_RCPP
}


SEXP bisq_knn_C( SEXP dd, SEXP hh ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXd d=Rcpp::as<VectorXd>(dd);
int h = Rcpp::as<int>(hh);
Function bisq_C("bisq_C");
d=d.cwiseAbs();
VectorXd rk = d;
std::sort(rk.data(),rk.data()+rk.size());
double dmax=rk(h);
while((rk(h+1)==dmax) | (dmax==0)){ h=h+1; dmax=rk(h)+0.00000001;}
SEXP aa=bisq_C(d,dmax,h);
VectorXd res = as<VectorXd>(aa);
return wrap(res);
END_RCPP
}


SEXP bin_C( SEXP dd, SEXP hh, SEXP Minv ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;
VectorXd dr=Rcpp::as<VectorXd>(dd);
VectorXd d=dr.array().abs();
double h = Rcpp::as<double>(hh);
int minv = Rcpp::as<int>(Minv);
Function sort("sort");
SEXP rkk =sort(d);
VectorXd zz=as<VectorXd>(rkk);
if(h<zz[minv]) h=zz[minv];
int n = d.size();
VectorXd z = d;
for(int j = 0; j < n; ++j)
{
if(d[j]>=h) {z[j] = 0;} else {z[j] = 1;}
}
return wrap(z);
END_RCPP
}


SEXP gauss_C( SEXP dd, SEXP hh, SEXP Minv ){
BEGIN_RCPP

 using namespace Rcpp;
 using Eigen::Map;
 using Eigen::VectorXd;
 using Eigen::Dynamic;

 VectorXd d = Rcpp::as<VectorXd>(dd);
 int minv = Rcpp::as<int>(Minv);
 double h = Rcpp::as<double>(hh);
 VectorXd y = d;
 std::sort(y.data(),y.data()+y.size());
 if (h < y(minv-1)) { h = y(minv-1);}

 double c =-0.5;
 VectorXd x = d/h;
 x=(x.array().square()*c).array().exp();
 return wrap(x);

END_RCPP
}


SEXP gauss_knn_C( SEXP dd, SEXP hh ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXd d=Rcpp::as<VectorXd>(dd);
int h = Rcpp::as<int>(hh);
Function EDK("EDK");
Function sort("sort");
bool check;
d=d.cwiseAbs();
SEXP rkk =sort(d);
VectorXd rk=as<VectorXd>(rkk);
double dmax=rk(h);
while((rk(h+1)==dmax) | (dmax==0)){ h=h+1;dmax=rk(h)+0.00000001;}
int n=d.size();
SEXP aa=EDK(d,0.5*dmax,2);
VectorXd res = as<VectorXd>(aa);
for(int j = 0; j < n; ++j)
{  check = (d(j) > dmax);
if(check) {res(j)=0;}
}
return wrap(res);
END_RCPP
}


SEXP gauss_adapt_C( SEXP dd, SEXP hh ){
BEGIN_RCPP

 using namespace Rcpp;
 using Eigen::Map;
 using Eigen::VectorXd;
 using Eigen::Dynamic;

 VectorXd d = Rcpp::as<VectorXd>(dd);
 int h = Rcpp::as<int>(hh);
 VectorXd y = d;
 std::sort(y.data(),y.data()+y.size());
 double dmax=y(h);
 double c =-0.5;
 VectorXd x = d/dmax;
 x=(x.array().square()*c).array().exp();
 return wrap(x);

END_RCPP
}


SEXP Dx_dense_C( SEXP xxi, SEXP xx, SEXP TIME ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
VectorXd x=Rcpp::as<VectorXd>(xx);
const double xi= Rcpp::as<double>(xxi) ;
bool time = as<bool>(TIME);
VectorXd x1= x.array()-xi;
VectorXd res= x1.array().square();
res = res.array().sqrt();
if(time) {return wrap(x1);} else {return wrap(res);}
END_RCPP
}


SEXP set0_conditionXD( SEXP XX, SEXP ZZ ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;

VectorXd x=Rcpp::as<VectorXd>(XX);
VectorXd z=Rcpp::as<VectorXd>(ZZ);

int n=x.size();
for(int j = 0; j < n; ++j)
{
if(z[j]>=0) {x[j] = 0;}
}
return wrap(x);
END_RCPP
}


SEXP kernelW_C( SEXP XX, SEXP hh, SEXP MykernelS, SEXP isgcv_, SEXP Type, SEXP Minv, SEXP maxknn_, SEXP NmaxDist_, SEXP TIME, SEXP Decay, SEXP DDiagNull ){
BEGIN_RCPP

using namespace Rcpp;
 using Eigen::Map;
 using Eigen::SparseMatrix;
 using Eigen::VectorXd;
  using Eigen::VectorXi;
 using Eigen::MatrixXd;
 using Eigen::MatrixXi;
 using Eigen::SparseVector;
 using Eigen::Dynamic;
 typedef Eigen::Triplet<double> T;

// input var
StringVector mykernels=as<StringVector>(MykernelS);
const String type=Rcpp::as<String>(Type);
MatrixXd X=Rcpp::as<MatrixXd>(XX);
VectorXd H = Rcpp::as<VectorXd>(hh);
double decay = Rcpp::as<double>(Decay);
int minv = Rcpp::as<int>(Minv);
int maxknn=Rcpp::as<int>(maxknn_);
int NmaxDist=Rcpp::as<int>(NmaxDist_);
bool isgcv = as<bool>(isgcv_);
bool time = as<bool>(TIME);
bool DiagNull = as<bool>(DDiagNull);
//bool Wnorm = as<bool>(normWW);
int n=X.rows();

//functions
Function GPKj("GPKj");


SparseVector<double> ds;
typedef Eigen::Triplet<double> T;
std::vector<T> tripletList;
SparseMatrix<double> W(n,n);
VectorXd z;
SEXP zz;
for(int j=0;j<n;++j){
	zz=GPKj(j,0,0,X,H,mykernels,type,minv,maxknn,NmaxDist,isgcv,time,decay);
	z = as<VectorXd>(zz);
	if(DiagNull) z(j)=0;
	ds = z.sparseView();
	for (SparseVector<double>::InnerIterator it(ds); it; ++it){
		tripletList.push_back(T(j,it.index(),it.value()));
	}
}
W.setFromTriplets(tripletList.begin(), tripletList.end());
return wrap(W) ;

END_RCPP
}


SEXP kernel_C( SEXP XX, SEXP J, SEXP hh, SEXP Mykernel, SEXP Minv, SEXP TIME, SEXP Decay, SEXP DDiagNull, SEXP normWW ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

const String mykernel=Rcpp::as<String>(Mykernel);
MatrixXd X=Rcpp::as<MatrixXd>(XX);
int j = Rcpp::as<int>(J);
double h = Rcpp::as<double>(hh);
double decay = Rcpp::as<double>(Decay);
int minv = Rcpp::as<int>(Minv);
bool time = as<bool>(TIME);
bool DiagNull = as<bool>(DDiagNull);
bool Wnorm = as<bool>(normWW);



//functions
Function bisq_knn_C("bisq_knn_C");
Function bin_C("bin_C");
Function rep("rep");
Function bisq_C("bisq_C");
Function gauss_C("gauss_C");
Function gauss_knn_C("gauss_knn_C");
Function gauss_adapt_C("gauss_adapt_C");
Function set0_conditionXD("set0_conditionXD");
Function fnormW("normW");
Function cbind("cbind");
Function D_dense_C("D_dense_C");
Function Dx_dense_C("Dx_dense_C");
SEXP  zz;
VectorXd z;
VectorXd d;
VectorXd res;

SEXP dd;

// if(X.cols()>2) {dd=maxdist_C(distances,j);}
if(X.cols()==2) {dd=D_dense_C(X(j,0),X(j,1),X.col(0),X.col(1));} else {dd=Dx_dense_C(X(j),X,time);}
d=as<VectorXd>(dd);

VectorXd delay=as<VectorXd>(rep(decay,X.rows()));
d=d+delay;

if(mykernel=="gauss")  {zz = gauss_C(d,h,minv);}  else if  (mykernel=="gauss_adapt")  {zz = gauss_adapt_C(d,h);}  else if  (mykernel=="gauss_knn")  {zz = gauss_knn_C(d,h);}  else if(mykernel=="bin") { zz = bin_C(d,h,minv);} else if(mykernel=="bisq_knn") { zz = bisq_knn_C(d,h);}
else if(mykernel=="bisq"){zz = bisq_C(d,h,minv);}

if(time){
SEXP z1=set0_conditionXD(zz,d);
z=as<VectorXd>(z1);
} else {z = as<VectorXd>(zz);}
if(DiagNull) {z(j)=0;}

if((Wnorm)&(z.sum()>0)) {res=z/z.sum();return wrap(res);} else {return wrap(z);}
//return wrap(dd);

END_RCPP
}


SEXP PhWY_C( SEXP YY, SEXP XX, SEXP WW, SEXP Wi ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

const Map<VectorXd> Y(as<Map<VectorXd> > (YY)); // Y peut être global ou local Y=Y*wi;
const Map<MatrixXd> X(as<Map<MatrixXd> >(XX)); // X peut être global ou local X=X*wi;
const SparseMatrix<double> W(as<SparseMatrix<double> >(WW));
const Map<VectorXd> wi(as<Map<VectorXd> > (Wi));
Function INST_C("INST_C");
Function Proj_C("Proj_C");

SEXP HH= INST_C(X,W,false,0);
MatrixXd H = as<MatrixXd>(HH);
H = H.array()*((wi.replicate(1,H.cols())).array());
VectorXd YY=Y.array()*wi.array();
MatrixXd WY=W*YY;
SEXP PhWY =Proj_C(H,WY);
MatrixXd PHWY=as<VectorXd>(PhWY);
return wrap(PHWY);

END_RCPP
}


SEXP fastlmLLT_C( SEXP XX, SEXP YY, SEXP SE_ ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;

MatrixXd X=Rcpp::as<MatrixXd>(XX);
MatrixXd Y=Rcpp::as<MatrixXd>(YY);
bool SE = as<bool>(SE_);
int n=Y.rows();
int p;
double s;
int df;
VectorXd fitted;
VectorXd resid;
VectorXd se;

VectorXd betahat;

const LLT<MatrixXd> llt(AtA(X));
 	betahat = llt.solve(X.adjoint() * Y);
 	//new
 	if(SE) {
    fitted=X * betahat;
    resid=Y - fitted;
    p=X.cols();
    df=n - p;
    s=resid.norm() /std::sqrt(double(df));
    se=s * llt.matrixL().solve(MatrixXd::Identity(p, p)).colwise().norm();
 	return List::create(Named("Betav")=betahat,Named("se")=se);
 	} else {
 	return List::create(Named("Betav")=betahat);
 	}

END_RCPP
}


SEXP QRcpp_C( SEXP AA, SEXP bb ){
BEGIN_RCPP

  using Eigen::Map;
  using Eigen::MatrixXd;
  typedef Eigen::ColPivHouseholderQR<MatrixXd> QR;
  MatrixXd A(as<MatrixXd>(AA));
  MatrixXd b(as<MatrixXd>(bb));
  MatrixXd X;
  QR solverQR(A);
  solverQR.compute(A);
  X = solverQR.solve(b);
  return wrap(X);
END_RCPP
}


SEXP QRcpp2_C( SEXP AA, SEXP bb, SEXP cc ){
BEGIN_RCPP

  using Eigen::Map;
  using Eigen::MatrixXd;
  typedef Eigen::ColPivHouseholderQR<MatrixXd> QR;
  MatrixXd A(as<MatrixXd>(AA));
  MatrixXd b(as<MatrixXd>(bb));
  MatrixXd c(as<MatrixXd>(cc));

  MatrixXd SY;
  MatrixXd XCw;

  QR solverQR(A);
  solverQR.compute(A);
  SY = solverQR.solve(b);
  XCw = solverQR.solve(c);

  return List::create(Named("SY")=SY,Named("XCw")=XCw);
END_RCPP
}


SEXP ApproxiW( SEXP WW, SEXP la, SEXP order ){
BEGIN_RCPP

using Eigen::SparseMatrix;
 using namespace Rcpp;
 using Eigen::MatrixXd;
 int ord = Rcpp::as<int>(order) ;
 double lambda = Rcpp::as<double>(la) ;
 const SparseMatrix<double> W(as<SparseMatrix<double> >(WW));
 int n = W.rows();
 SparseMatrix<double> A = W;
 double b = lambda;
 typedef Eigen::Triplet<double> T;
 std::vector<T> tripletList;
 tripletList.reserve(n);
 for(int j = 0; j < n; ++j)
{
 tripletList.push_back(T(j,j,1));
}
 SparseMatrix<double> iW(n,n);
 iW.setFromTriplets(tripletList.begin(), tripletList.end());
 iW = iW+lambda*W;
 for(int j = 2; j < ord; ++j)
{
 A = A * W;
 b = lambda*b;
 iW = iW + b*A;
}
 //SparseMatrix<double> iWW = iW.transpose()*iW;
 //return List::create(Named("iW")=iW,Named("iWW")=iWW);
 return wrap(iW);
END_RCPP
}


SEXP mod( SEXP YY, SEXP XX, SEXP WW, SEXP XZZ, SEXP YZZ, SEXP Wi, SEXP LocalInst, SEXP ismethodB2SLS, SEXP ismethodMGWRSAR_1_kc_0, SEXP SE_ ){
BEGIN_RCPP

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

const Map<VectorXd> Y(as<Map<VectorXd> > (YY)); // Y peut être global ou local Y=Y*wi;
const Map<MatrixXd> X(as<Map<MatrixXd> >(XX)); // X peut être global ou local X=X*wi;
const SparseMatrix<double> W(as<SparseMatrix<double> >(WW));
const Map<VectorXd> YZ(as<Map<VectorXd> > (YZZ)); // YZ=Y est global
const Map<MatrixXd> XZ(as<Map<MatrixXd> >(XZZ)); // XZ=XZ est global
const Map<VectorXd> wi(as<Map<VectorXd> > (Wi));
const String localInst=Rcpp::as<String>(LocalInst);
bool isB2SLS = as<bool>(ismethodB2SLS);
bool isMGWRSAR_1_kc_0 = as<bool>(ismethodMGWRSAR_1_kc_0);
bool SE = as<bool>(SE_);

Function INST_C("INST_C");
Function Proj_C("Proj_C");
Environment Mat("package:Matrix");
Function cBind = Mat["cBind"];

MatrixXd WY;
SEXP HH;
MatrixXd H;
SEXP PhWY;
MatrixXd XB;
SEXP XXB;
double lambda;
VectorXd l;
VectorXd betahat;
VectorXd res;

////////////////////////////////////////////////////////////////
// si localInst=L0 on utilise PhWY = [Z inv(Z tZ) tZ WY]
// si localInst=L1 on utilise PhWY = wi * [Z inv(Z tZ) tZ WY]
// si localInst=L2 on utilise PhWY = Z inv(Z tZ) tZ W wi*Y
// avec Zi = [ X*wi  W  X*wi  W^2 X*wi  ...]
// si localInst=L3 on utilise PhWY = Zi inv(Zi tZi) tZi W wi*Y  // OPTIMAL pour typ_loc=poly
// si localInst=L4 on utilise PhWY = wi * [Zi inv(Zi tZi) tZi WY]
// avec Zi = wi * Z
// si localInst=L5 on utilise PhWY = Zi inv(Zi tZi) tZi W wi*Y   //SECOND BEST pour typ_loc= random et typ_loc=poly
// si localInst=L6 on utilise PhWY = wi * [Zi inv(Zi tZi) tZi WY] // OPTIMAL pour typ_loc= random
// si localInst=L7 on utilise PhWY = wi * [ZZi inv(ZZi tZZi) tZZi WY] avec ZZi = [XVi XC]
////////////////////////////////////////////////////////////////


if(localInst=="L0"){HH= INST_C(XZ,W,false,0);H = as<MatrixXd>(HH);WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L1"){HH= INST_C(XZ,W,false,0);H = as<MatrixXd>(HH);WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if(localInst=="L2"){HH= INST_C(XZ,W,false,0);H = as<MatrixXd>(HH);
VectorXd Yz=YZ.array()*wi.array();WY=W*Yz;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L3"){
MatrixXd Xz=XZ.array()*((wi.replicate(1,XZ.cols())).array());
HH= INST_C(Xz,W,false,0);H = as<MatrixXd>(HH);
VectorXd Yz=YZ.array()*wi.array();WY=W*Yz;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L4"){
MatrixXd Xz=XZ.array()*((wi.replicate(1,XZ.cols())).array());
HH= INST_C(Xz,W,false,0);H = as<MatrixXd>(HH);
WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if(localInst=="L5"){
HH= INST_C(XZ,W,false,0);
H = as<MatrixXd>(HH);
H = H.array()*((wi.replicate(1,H.cols())).array());
VectorXd Yz=YZ.array()*wi.array();
WY=W*Yz;
PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L6"){
HH= INST_C(XZ,W,false,0);H = as<MatrixXd>(HH);H = H.array()*((wi.replicate(1,H.cols())).array());
WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if(localInst=="L7"){
HH= INST_C(XZ,W,false,0);H = as<MatrixXd>(HH);
SEXP wi1=cBind(XZ.col(1),XZ.col(1),wi,wi);
MatrixXd wi2=as<MatrixXd>(wi1);
H = H.array()*(wi2.array());
WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}


if (!isMGWRSAR_1_kc_0) {XXB = cBind(X,XB);XB= as<MatrixXd>(XXB);}

const LLT<MatrixXd> llt(AtA(XB));
betahat = llt.solve(XB.adjoint() * Y);
l = betahat.tail(1);
lambda = l(0);



if(isB2SLS) {
if(localInst=="L0"){HH= INST_C(XZ,W,true,lambda);H = as<MatrixXd>(HH);WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L1"){HH= INST_C(XZ,W,true,lambda);H = as<MatrixXd>(HH);WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if(localInst=="L2"){HH= INST_C(XZ,W,true,lambda);H = as<MatrixXd>(HH);
VectorXd Yz=YZ.array()*wi.array();WY=W*Yz;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L3"){
MatrixXd Xz=XZ.array()*((wi.replicate(1,XZ.cols())).array());
HH= INST_C(Xz,W,true,lambda);H = as<MatrixXd>(HH);
VectorXd Yz=YZ.array()*wi.array();WY=W*Yz;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L4"){
MatrixXd Xz=XZ.array()*((wi.replicate(1,XZ.cols())).array());
HH= INST_C(Xz,W,true,lambda);H = as<MatrixXd>(HH);
WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if(localInst=="L5"){
HH= INST_C(XZ,W,true,lambda);H = as<MatrixXd>(HH);H = H.array()*((wi.replicate(1,H.cols())).array());VectorXd Yz=YZ.array()*wi.array();WY=W*Yz;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY);}

if(localInst=="L6"){
HH= INST_C(XZ,W,true,lambda);H = as<MatrixXd>(HH);H = H.array()*((wi.replicate(1,H.cols())).array());
WY=W*YZ;PhWY =Proj_C(H,WY);
XB=as<VectorXd>(PhWY).array()*wi.array();}

if (!isMGWRSAR_1_kc_0) {XXB = cBind(X,XB);XB= as<MatrixXd>(XXB);}
const LLT<MatrixXd> llt2(AtA(XB));
betahat = llt2.solve(XB.adjoint() * Y);
}

int n=Y.rows();
int p;
double s;
int df;
VectorXd fitted;
VectorXd resid;
VectorXd se;

if(SE) {
    fitted=X * betahat;
    resid=Y - fitted;
    p=X.cols();
    df=n - p;
    s=resid.norm() /std::sqrt(double(df));
    se=s * llt.matrixL().solve(MatrixXd::Identity(p, p)).colwise().norm();
 	return List::create(Named("Betav")=betahat,Named("se")=se);
 	} else {
 	return List::create(Named("Betav")=betahat);
 	}

END_RCPP
}
