// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec cpp_cr_rho(arma::vec x, double theta) {
  arma::vec xv(x.size());

  if(theta == 0){
    xv = -exp(-x);
  }else if(theta == -1){
    xv = log(1+x);
  }else{
    xv = -1*pow((1-theta*x), 1+1/theta);
    xv = xv/(1+theta);
  }

  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);

  return xv;
}


// [[Rcpp::export]]
arma::vec cpp_dcr_rho(arma::vec x, double theta) {
  arma::vec xv(x.size());


  if(theta == 0){
    xv = exp(-x);
  }else if(theta == -1){
    xv = 1/(1+x);
  }else{
    xv = pow((1-theta*x),1/theta);
  }

  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);

  return xv;
}


// [[Rcpp::export]]
arma::vec cpp_ddcr_rho(arma::vec x, double theta) {
  arma::vec xv(x.size());

  //double intpart;

  //return modf(theta, &intpart) == 0.0;

  if(theta == 0){
    xv = -exp(-x);
  }else if(theta == -1){
    xv = -1/pow((1+x),2);
  }else{
    xv = -1*pow((1-theta*x),1/theta-1);
  }

  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);

  return xv;
}


// [[Rcpp::export]]
double cpp_obj(NumericVector lam, NumericMatrix u, NumericVector ubar,
               arma::vec Ti, double theta){

  arma::mat umat(u.begin(), u.nrow() , u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  arma::vec lam_t_u = umat*lambda;
  double lam_t_ubar = as_scalar(u_bar*lambda);

  return -mean(na_omit(as<NumericVector>(wrap(Ti % cpp_cr_rho( lam_t_u, theta))))) + lam_t_ubar;
}


// [[Rcpp::export]]
arma::rowvec cpp_derv_obj(NumericVector lam,NumericMatrix u, NumericVector ubar,
                    arma::vec Ti, double theta){
  int N = u.nrow();
  arma::mat umat(u.begin(), N , u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  arma::vec lam_t_u = umat*lambda;
  arma::rowvec temp = trans(Ti % cpp_dcr_rho( lam_t_u, theta));
  temp.elem(find_nonfinite(temp)).fill(0);


  return u_bar-( temp*umat)/N;
}


// [[Rcpp::export]]
arma::mat cpp_update_hessianInv(arma::mat oldInv, arma::vec sk, arma::vec yk){
  double scal1 = as_scalar(trans(yk)*oldInv*yk);
  double scal2 = as_scalar(trans(sk)*yk);
  arma::mat mat1 = (oldInv*yk)*trans(sk);

  return oldInv + ((scal1+scal2)/pow(scal2,2))*sk*trans(sk) - (mat1 + trans(mat1))/scal2;
}


// [[Rcpp::export]]
double cpp_backtrack(double alpha,double beta, NumericVector x_val,
                     NumericVector del_x,NumericVector nabla_f,
                     NumericMatrix u, NumericVector ubar, arma::vec Ti,
                       double theta){

  int N = u.nrow();
  arma::rowvec nablaf(nabla_f.begin(), nabla_f.size(), false);
  arma::colvec delx(del_x.begin(), del_x.size(), false);

  arma::mat umat(u.begin(), N , u.ncol(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);


  double step = 1.0;

  double f_x = cpp_obj(x_val, u, ubar, Ti,theta);
  double df_t_dx =  as_scalar(nablaf*delx);

  while(cpp_obj(x_val+step*del_x,u,ubar,Ti,theta) > f_x+alpha*step*df_t_dx ){
    step = beta*step;
  }
  return step;

}

// [[Rcpp::export]]
List cpp_quasi_newt(NumericVector ini,NumericMatrix u, NumericVector ubar,
                    arma::vec Ti, double theta, double alpha, double beta,
                    int max_iter, double tol){
  int p = ini.size();
  int N = u.nrow();
  arma::mat umat(u.begin(), N , u.ncol(), false);

  arma::vec current_est(ini.begin(), ini.size(), false);
  arma::vec new_est;


  arma::vec current_derv, new_derv;
  arma::vec current_direction, new_direction;

  arma::mat current_InvHessian = eye<mat>(p,p);
  arma::mat new_InvHessian;

  double stepSize, objValue;

  arma::vec yk, sk;

  for(int i=1; i<= max_iter; i++){
    current_derv = trans( cpp_derv_obj(as<NumericVector>(wrap(current_est)), u,ubar,Ti,theta) );
    current_direction = -current_InvHessian*current_derv;

    objValue = cpp_obj(as<NumericVector>(wrap(current_est)), u,  ubar, Ti,  theta);
    //cout<< objValue << "\n";
    stepSize = cpp_backtrack(alpha,beta, as<NumericVector>(wrap(current_est)),
                  as<NumericVector>(wrap(current_direction)),
                  as<NumericVector>(wrap(current_derv)),
                  u,ubar, Ti,theta);
    //cout<< stepSize << "\n";
    new_est = current_est + stepSize*current_direction;
    //cout<< new_est;

    new_derv = trans(cpp_derv_obj(as<NumericVector>(wrap(new_est)), u,ubar,Ti,theta));
    yk = new_derv - current_derv;
    sk = new_est - current_est;
    new_InvHessian = cpp_update_hessianInv(current_InvHessian, sk, yk);


    if(objValue< -1e+30){
      Function warning("warning");
      warning("The objective function is unbounded, a different theta might be needed.");
      arma::vec weights = zeros<vec>(N);
      arma::uvec indx = find(Ti == 1);
      arma::vec lambda(new_est.begin(), new_est.size(), false);
      weights.elem(indx) = (cpp_dcr_rho( umat*lambda, theta)).elem(indx)/N;
      return List::create(Named("res") = new_est, Named("weights") = weights,
                          Named("Conv") = false);
    }

    //cout<< sum(square(new_derv))<< "\n";
    if(sum(square(new_derv)) < tol ){
      arma::vec weights = zeros<vec>(N);
      arma::uvec indx = find(Ti == 1);
      arma::vec lambda(new_est.begin(), new_est.size(), false);
      weights.elem(indx) = (cpp_dcr_rho( umat*lambda, theta)).elem(indx)/N;
      return List::create(Named("res") = new_est, Named("weights") = weights,
                          Named("Conv") = true);
    }else{
      current_est = new_est;
      current_InvHessian = new_InvHessian;

    }



  }
  arma::vec weights = zeros<vec>(N);
  arma::uvec indx = find(Ti == 1);
  arma::vec lambda(new_est.begin(), new_est.size(), false);
  weights.elem(indx) = (cpp_dcr_rho( umat*lambda, theta)).elem(indx)/N;

  return List::create(Named("res") = new_est, Named("weights") = weights,
                      Named("Conv") = false);
}




// // [[Rcpp::export]]
// arma::mat get_gk_simple(NumericMatrix u, arma::vec Y, arma::vec Ti,
//                         arma::vec lam , arma::vec beta,
//                         double tau1, double tau0, double theta){
//   int N = Y.size();
//   arma::mat umat(u.begin(), N , u.ncol(), false);
//   arma::mat ucopy1(u.begin(), N , u.ncol(), true);
//   arma::mat ucopy2(u.begin(), N , u.ncol(), true);
//
//   arma::vec temp1 = Ti % cpp_dcr_rho(umat*lam, theta);
//   ucopy1.each_col() %= temp1;
//   arma::mat gk1 = ucopy1 - umat;
//
//   arma::vec temp2 = (1-Ti) % cpp_dcr_rho(umat*beta, theta);
//   ucopy2.each_col() %= temp2;
//   arma::mat gk2 = ucopy2 - umat;
//
//   arma::vec gk3 = (Ti % Y)%cpp_dcr_rho(umat*lam,theta) - tau1;
//   arma::vec gk4 = ((1-Ti) % Y)%cpp_dcr_rho(umat*beta,theta) - tau0;
//   return join_rows(join_rows(join_rows(gk1,gk2),gk3),gk4);
// }
//



// // [[Rcpp::export]]
// arma::mat cpp_derv2_obj(NumericVector lam,NumericMatrix u, NumericMatrix ucopy, NumericVector ubar,
//                           arma::vec Ti, double theta){
//   int N = u.nrow();
//   arma::mat umat(u.begin(), N , u.ncol(), false);
//   //I create a copy here which I need to use later.
//   //uses extra memory but I don't see a way around it yet.
//   arma::mat umatCopy(ucopy.begin(), N , ucopy.ncol(), false);
//
//   arma::colvec lambda(lam.begin(), lam.size(), false);
//   arma::rowvec u_bar(ubar.begin(), ubar.size(), false);
//
//   arma::vec lam_t_u = umat*lambda;
//   arma::rowvec temp = trans(Ti % cpp_ddcr_rho( lam_t_u, theta));
//   temp.elem(find_nonfinite(temp)).fill(0);
//
//   umatCopy.each_col() %= trans(temp);
//   return -(trans(umat)*umatCopy)/N;
// }
