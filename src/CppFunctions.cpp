// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

//A function for evaluating the Cressie Read family
//See documentation for the exact form of functions
//In package Vignette this function is used for the function \rho
//Note: In this package we use the functions with 1-theta as opposed
// to 1+theta.
// [[Rcpp::export]]
arma::vec cpp_cr_rho(arma::vec x, double theta) {
  //Create a copy of the vector
  arma::vec xv(x.size());


  if(theta == 0){
    xv = -exp(x); //The limiting case of theta = 0. Easy to calculate
  }else if(theta == -1){
    xv = log(1-x);// The limiting case of theta = -1, use L'Hopital's rule
  }else{
    // The function we use for all other values of theta
    xv = -1*pow((1+theta*x), 1+1/theta);
    xv = xv/(1+theta);
  }

  //For some cases particularly theta == -1, the function is
  //not defined for the entire real line. In order to maintain
  //concavity of the function we assign the value -Inf to regions
  //outside the domain of the function.
  //This is a standard trick in convex optimization
  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);
  return xv;
}

//The first derivative of the above function evaluated for a vector x
// [[Rcpp::export]]
arma::vec cpp_dcr_rho(arma::vec x, double theta) {
  arma::vec xv(x.size());

  if(theta == 0){
    xv = -exp(x);
  }else if(theta == -1){
    xv = -1/(1-x);
  }else{
    xv = -1*pow((1+theta*x),1/theta);
  }

  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);

  return xv;
}

//The second derivative of the above function evaluated for a vector x
// [[Rcpp::export]]
arma::vec cpp_ddcr_rho(arma::vec x, double theta) {
  arma::vec xv(x.size());

  if(theta == 0){
    xv = -exp(x);
  }else if(theta == -1){
    xv = -1/pow((1-x),2);
  }else{
    xv = -1*pow((1+theta*x),1/theta-1);
  }

  xv.elem(find_nonfinite(xv)).fill(-1*datum::inf);

  return xv;
}

//////////////////////////////////////////////////////////

//The main objective function which we need to optimize over.
//See package Vignette for the objective functions/optimization problem
//
//Note that R is pretty flexible with making the conversions between R objects
// and Rcpp objects. E.g. a vector in R can be used for the following inputs:
// NumericVector arma::vec arma::colvec arma::rowvec
// [[Rcpp::export]]
double cpp_obj(NumericVector lam, NumericMatrix u, NumericVector ubar,
               arma::vec Ti, double theta){

  //Convert some Rcpp objects into armadillo objects
  //The false allows us to do this quickly and withou alloting
  //extra memory. HOWEVER, since this just assigns a pointer
  //making changes to umat will also change u
  arma::mat umat(u.begin(), u.nrow() , u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  //This creates a vector where the i-th term is
  // \lambda^T u_K(X_i)
  arma::vec lam_t_u = umat*lambda;
  double lam_t_ubar = as_scalar(u_bar*lambda);

  // Return the value of the objective function.
  //Omit any NA terms caused by Ti==0 and cpp_cr_rho( lam_t_u, theta) == -Inf
  //as<>wrap() is needed to use the na_omit sugar function since it is an Rcpp
  //function and not available in RcppAramdillo
  // The operator % is the element wise /Schur product
  return -mean(na_omit(as<NumericVector>(wrap(Ti % cpp_cr_rho( lam_t_u, theta))))) + lam_t_ubar;
}


//The first derivative of the above objective function
//This function returns a row vector of size K where K = ncol(u)
// [[Rcpp::export]]
arma::rowvec cpp_derv_obj(NumericVector lam,NumericMatrix u, NumericVector ubar,
                    arma::vec Ti, double theta){
  //Initialize the vectors
  //row vs. col vec is declared for computational ease only
  //to allow for proper matrix operations
  int N = u.nrow();
  arma::mat umat(u.begin(), N , u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  //Same way as before we create a vector for the terms
  // \lambda^T u_K(X_i)
  arma::vec lam_t_u = umat*lambda;
  arma::rowvec temp = trans(Ti % cpp_dcr_rho( lam_t_u, theta));
  temp.elem(find_nonfinite(temp)).fill(0);

  //Return derivative
  return u_bar-( temp*umat)/N;
}


//This function gives us a one-step update for the BFGS algorithm
//oldInv is the current "Approximate Inverse Hessian"
//sk is the difference x_{k+1} - x_{k}
//yk is the difference f'(x_{k+1}) - f'(x_{k})
//where x_k is the optimization variable and f is the objective function
// [[Rcpp::export]]
arma::mat cpp_update_hessianInv(arma::mat oldInv, arma::vec sk, arma::vec yk){
  //We evaluate some scalars first since we use them in multiple locations
  // of the updating step
  double scal1 = as_scalar(trans(yk)*oldInv*yk);
  double scal2 = as_scalar(trans(sk)*yk);
  arma::mat mat1 = (oldInv*yk)*trans(sk);

  //The BFGS update step here:
  //Formulation from https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm
  return oldInv + ((scal1+scal2)/pow(scal2,2))*sk*trans(sk) - (mat1 + trans(mat1))/scal2;
}


//A simple backtracking lin search
//This formulation/notation is taken from Algorithm 9.2 of
//Boyd, Stephen, and Lieven Vandenberghe. Convex optimization. 2004.

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

  //Begin with an initial step size of 1
  double step = 1.0;

  //Initialize some quantities which do not change
  double f_x = cpp_obj(x_val, u, ubar, Ti,theta);
  double df_t_dx =  as_scalar(nablaf*delx);

  //main loop for the algoritm to find an appropriate step size
  //Loop until the objective function of the new proposed step is sufficiently small
  while(cpp_obj(x_val+step*del_x,u,ubar,Ti,theta) > f_x+alpha*step*df_t_dx ){
    step = beta*step;
  }
  return step;
}

//The main function for obtaining the point estimates.
//This function implements a BFGS algorithm

// [[Rcpp::export]]
List cpp_quasi_newt(NumericVector ini,NumericMatrix u, NumericVector ubar,
                    arma::vec Ti, double theta, double alpha, double beta,
                    int max_iter, double tol){

  //Initialize some variables
  int p = ini.size();
  int N = u.nrow();
  arma::mat umat(u.begin(), N , u.ncol(), false);

  //Initialize variables which will be used in the main loop
  arma::vec current_est(ini.begin(), ini.size(), false);
  arma::vec new_est;

  arma::vec current_derv, new_derv;
  arma::vec current_direction, new_direction;

  arma::mat current_InvHessian = eye<mat>(p,p);
  arma::mat new_InvHessian;

  double stepSize, objValue;

  arma::vec yk, sk;

  for(int i=1; i<= max_iter; i++){
    //Calculate the current derivative of the objective function and direction
    current_derv = trans( cpp_derv_obj(as<NumericVector>(wrap(current_est)), u,ubar,Ti,theta) );
    current_direction = -current_InvHessian*current_derv;

    //Evaluate the objective function as the current point
    //We use this later to check if the objective is unbounded
    objValue = cpp_obj(as<NumericVector>(wrap(current_est)), u,  ubar, Ti,  theta);

    //Calculate the Stepsize using the backtracking line search algorithm
    stepSize = cpp_backtrack(alpha,beta, as<NumericVector>(wrap(current_est)),
                  as<NumericVector>(wrap(current_direction)),
                  as<NumericVector>(wrap(current_derv)),
                  u,ubar, Ti,theta);

    //Calulate the new estimate
    new_est = current_est + stepSize*current_direction;

    //Update the estimate for the inverse hessian using the BFGS update
    new_derv = trans(cpp_derv_obj(as<NumericVector>(wrap(new_est)), u,ubar,Ti,theta));
    yk = new_derv - current_derv;
    sk = new_est - current_est;
    new_InvHessian = cpp_update_hessianInv(current_InvHessian, sk, yk);

    //For some cases, say theta = -1, the objective function can be unbounded
    //in which case the alorithm will never converge.
    //In this case we throw a warning suggesting using a different theta value
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

    //Check convergence by checking the first derivative
    //Since the derivative would be 0 at a critical point we
    //calculate the l2 norm to check convergence
    if(sum(square(new_derv)) < tol ){
      //If the stopping condition is met
      //Calculate the weights used for getting our final estimate
      arma::vec weights = zeros<vec>(N);
      arma::uvec indx = find(Ti == 1);
      arma::vec lambda(new_est.begin(), new_est.size(), false);
      weights.elem(indx) = (cpp_dcr_rho( umat*lambda, theta)).elem(indx)/N;
      return List::create(Named("res") = new_est, Named("weights") = weights,
                          Named("Conv") = true);
    }else{
      //If stopping condition is not met update the current estimate and
      //current InvHessian
      current_est = new_est;
      current_InvHessian = new_InvHessian;
    }
  }

  //If algotihm runs for maximum number of iterations
  //return the last result and state convergence status is false
  arma::vec weights = zeros<vec>(N);
  arma::uvec indx = find(Ti == 1);
  arma::vec lambda(new_est.begin(), new_est.size(), false);
  weights.elem(indx) = (cpp_dcr_rho( umat*lambda, theta)).elem(indx)/N;

  return List::create(Named("res") = new_est, Named("weights") = weights,
                      Named("Conv") = false);
}
