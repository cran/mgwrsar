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
  SEXP Sl_C( SEXP llambda, SEXP WW, SEXP iinv, SEXP aapprox) ;
  SEXP INST_C( SEXP XX, SEXP WW, SEXP withlambda, SEXP llambda) ;
  //SEXP set0_conditionXD( SEXP XX, SEXP ZZ) ;
  SEXP PhWY_C( SEXP YY, SEXP XX, SEXP WW, SEXP Wi) ;
  SEXP QRcpp2_C( SEXP AA, SEXP bb, SEXP cc) ;
  SEXP ApproxiW( SEXP WW, SEXP la, SEXP order) ;
  SEXP ApproxiWv( SEXP WW, SEXP la, SEXP order) ;
  SEXP mod( SEXP YY, SEXP XX, SEXP WW, SEXP XZZ, SEXP YZZ, SEXP Wi, SEXP LocalInst, SEXP ismethodB2SLS, SEXP ismethodMGWRSAR_1_kc_0, SEXP SE_) ;
}

static const R_CallMethodDef CallEntries[] = {
  {"ApproxiW",         (DL_FUNC) &ApproxiW,          3},
  {"ApproxiWv",        (DL_FUNC) &ApproxiWv,         3},
  {"QRcpp2_C",         (DL_FUNC) &QRcpp2_C,          3},
  {"INST_C",           (DL_FUNC) &INST_C,            4},
  {"mod",              (DL_FUNC) &mod,              10},
  {"PhWY_C",           (DL_FUNC) &PhWY_C,            4},
  {"Proj_C",           (DL_FUNC) &Proj_C,            2},
  //{"set0_conditionXD", (DL_FUNC) &set0_conditionXD,  2},
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
//Environment Mat("package:Matrix");
  //Function cBind = Mat["cBind"];

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
      H = cbind(X,WX,WWX,WWWX);
    }
    HH=as<MatrixXd>(H);
    return wrap(HH);
    END_RCPP
}


// SEXP set0_conditionXD( SEXP XX, SEXP ZZ ){
//   BEGIN_RCPP
//
//   using namespace Rcpp;
//   using Eigen::Map;
//   using Eigen::VectorXd;
//
//   VectorXd x=Rcpp::as<VectorXd>(XX);
//   VectorXd z=Rcpp::as<VectorXd>(ZZ);
//
//   int n=x.size();
//   for(int j = 0; j < n; ++j)
//   {
//     if(z[j]>=0) {x[j] = 0;}
//   }
//   return wrap(x);
//   END_RCPP
// }


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
  SparseMatrix<double> iW(n,n);
  iW.setIdentity();
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


SEXP ApproxiWv( SEXP WW, SEXP la, SEXP order ){
  BEGIN_RCPP

  using Eigen::SparseMatrix;
  using namespace Rcpp;
  using Eigen::MatrixXd;
  using namespace Eigen;
  using Eigen::VectorXd;
  Function row_prod("row_prod");

  int ord = Rcpp::as<int>(order) ;
  VectorXd lambda = Rcpp::as<VectorXd>(la);
  SparseMatrix<double> W = Rcpp::as<SparseMatrix<double> >(WW);
  int n = W.rows();
  SEXP BB= row_prod(lambda,W);
  SparseMatrix<double> B = as<SparseMatrix<double> >(BB);

  VectorXd b = lambda;
  SparseMatrix<double> iW(n,n);
  iW.setIdentity();
  iW = iW+B;
  SparseMatrix<double> A = W;
  for(int j = 2; j < ord; ++j)
  {
    A = A * W;
    b = lambda.array()*b.array();
    SEXP BB = row_prod(b,A);
    SparseMatrix<double> B = as<SparseMatrix<double> >(BB);
    iW = iW + B;
  }


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
  Function cbind("cbind");
  //Function cBind = Mat["cBind"];

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
    SEXP wi1=cbind(XZ.col(1),XZ.col(1),wi,wi);
    MatrixXd wi2=as<MatrixXd>(wi1);
    H = H.array()*(wi2.array());
    WY=W*YZ;PhWY =Proj_C(H,WY);
    XB=as<VectorXd>(PhWY).array()*wi.array();}


  if (!isMGWRSAR_1_kc_0) {XXB = cbind(X,XB);XB= as<MatrixXd>(XXB);}

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

    if (!isMGWRSAR_1_kc_0) {XXB = cbind(X,XB);XB= as<MatrixXd>(XXB);}
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
