// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec CRFamily(arma::vec x, double theta) {
  // A function for evaluating the Cressie Read family.
  // See documentation for the exact form for the class of functions.
  // Note: In this package we use a modification of this family
  //   Instead of 1 + theta we use 1 - theta.
  // Args:
  //   x: Vector of points at which function will be evaluated.
  //   theta: Scalar theta which parametrizes the CR family of
  //          functions.
  //
  // Returns:
  //   The vector with the CR function evaluated for each point
  //   in x.

  //Create a copy of the vector
  arma::vec xv(x.size());

  if (theta == 0) {
    xv = -exp(-x); //The limiting case of theta = 0. Easy to calculate.
  } else if(theta == -1) {
    xv = log(1 + x); // The limiting case of theta = -1, use L'Hopital's rule.
  } else {
    // The function we use for all other values of theta.
    xv = -1 * pow((1 - theta * x), 1 + 1 / theta);
    xv = xv / ( 1 + theta);
  }

  // For some cases particularly theta == -1, the function is
  // not defined for the entire real line. In order to maintain
  // concavity of the function we assign the value -Inf to regions
  // outside the domain of the function.
  // This is a standard trick in convex optimization.
  xv.elem(find_nonfinite(xv)).fill(-1 * datum::inf);
  return xv;
}

// [[Rcpp::export]]
arma::vec CRFamilyDerivative(arma::vec x, double theta) {
  // A function for evaluating the first derivative of
  // Cressie Read family.
  //
  // Args:
  //   x: Vector of points at which derivative will be evaluated.
  //   theta: Scalar theta which parametrizes the CR family of
  //          functions.
  //
  // Returns:
  //   The vector with the CR function's derivative evaluated
  //   for each point in x.

  arma::vec xv(x.size());

  if (theta == 0) {
    xv = exp(-x);
  } else if (theta == -1) {
    xv = 1 / (1 + x);
  } else {
    xv = pow(1 - theta * x, 1 / theta);
  }
  // Replace NA with -Inf as before
  xv.elem(find_nonfinite(xv)).fill(-1 * datum::inf);

  return xv;
}

// [[Rcpp::export]]
arma::vec CRFamilySecondDerivative(arma::vec x, double theta) {
  // A function for evaluating the second derivative of
  // the Cressie Read family.
  //
  // Args:
  //   x: Vector of points at which the second derivative will
  //      be evaluated.
  //   theta: Scalar theta which parametrizes the CR family of
  //          functions.
  //
  // Returns:
  //   The vector with the CR function's 2nd derivative evaluated
  //   for each point in x.

  arma::vec xv(x.size());

  // Consider the 3 cases again
  if (theta == 0) {
    xv = -exp(-x);
  } else if (theta == -1) {
    xv = -1 / pow(1 + x, 2);
  } else {
    xv = -1 * pow(1 - theta * x, 1 / theta - 1);
  }

  xv.elem(find_nonfinite(xv)).fill(-1 * datum::inf);

  return xv;
}


// [[Rcpp::export]]
double ObjectiveFunction(NumericVector lam,
                         NumericMatrix u,
                         NumericVector ubar,
                         arma::vec treat, double theta) {

  // The main objective function which we need to optimize over.
  // See package Vignette for the objective functions/optimization problem.
  //
  // Args:
  //   lam: Lambda vector at which to evaluate objective.
  //   u: A N * K matrix for the covariates. In our case this is
  //      this is just the matrix cbind(1, X) for a design matrix X.
  //      N is the number of subjects with K-1 covariates.
  //   ubar: The K vector of column means of matrix u.
  //   treat: The N vector specifying treatment assignment.
  //   theta: Scalar parametrizing the CR family of functions.
  //
  // Returns:
  //   A scalar value of the objective function evaluated at the vector
  //   lambda.


  // Convert some Rcpp objects into armadillo objects
  // The 'false' allows us to do this quickly and without alloting
  // extra memory. HOWEVER, since this just assigns a pointer
  // making changes to u_mat will also change u
  arma::mat u_mat(u.begin(), u.nrow(), u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  // This creates a vector where the i-th term is
  // \lambda^T * u_K(X_i)
  arma::vec lam_t_u = u_mat * lambda;
  double lam_t_ubar = as_scalar(u_bar * lambda);

  // Return the value of the objective function.
  // Omit any NA terms caused by treat == 0 and CRFamily(lam_t_u, theta) == -Inf
  // as<>wrap() is needed to use the na_omit sugar function since it is an Rcpp
  // function and not available in RcppAramdillo.
  // The operator % is the element wise/Schur product.
  return -mean(na_omit(as<NumericVector>(wrap(treat % CRFamily(lam_t_u, theta)))))
    + lam_t_ubar;
}


// [[Rcpp::export]]
arma::rowvec ObjectiveFirstDerivative(NumericVector lam,
                                      NumericMatrix u,
                                      NumericVector ubar,
                                      arma::vec treat, double theta) {

  // The first derivative of the main objective function
  // which we are minimizing here.
  //
  // Args:
  //   lam: Lambda vector at which to evaluate derivative.
  //   u: See ObjectiveFunction above.
  //   ubar: See ObjectiveFunction above.
  //   treat: See ObjectiveFunction above.
  //   theta: See ObjectiveFunction above.
  //
  // Returns:
  //   A scalar value of the objective function evaluated at the vector
  //   lambda.

  //Initialize the vectors.
  //Row vs. col vec is declared for computational ease only
  //to allow for proper matrix operations.
  int N = u.nrow();
  arma::mat u_mat(u.begin(), N, u.ncol(), false);
  arma::colvec lambda(lam.begin(), lam.size(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  //Same way as before we create a vector for the terms
  // \lambda^T u_K(X_i)
  arma::vec lam_t_u = u_mat * lambda;
  arma::rowvec temp = trans(treat % CRFamilyDerivative(lam_t_u, theta));
  temp.elem(find_nonfinite(temp)).fill(0);

  //Return derivative
  return u_bar - (temp * u_mat) / N;
}


// [[Rcpp::export]]
arma::mat UpdateInverseHessian(arma::mat old_inv_hessian,
                               arma::vec diff_in_est,
                               arma::vec diff_in_derv) {
  // This function gives us the one-step update for the BFGS algorithm
  // for the inverse of the Hessian matrix.
  //
  // Args:
  //   old_inv_hessian: The current "Approximate Inverse Hessian".
  //   diff_in_est: The difference between current and old estimate,
  //                x_{k+1} - x_{k}.
  //   diff_in_derv: The difference f'(x_{k+1}) - f'(x_{k})
  //                 where f is the objective function.
  // Returns:
  //   A matrix of the same dimensions as old_inv_hessian representing
  //   the one-step update for the inverse Hessian.

  // We evaluate some quantities which are used multiple times in the
  // computation.
  double scal1 = as_scalar(trans(diff_in_derv) * old_inv_hessian * diff_in_derv);
  double scal2 = as_scalar(trans(diff_in_est) * diff_in_derv);
  arma::mat mat1 = (old_inv_hessian * diff_in_derv) * trans(diff_in_est);

  //The BFGS update step :
  //Formulation from
  //https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm
  return old_inv_hessian + ((scal1 + scal2) / pow(scal2, 2)) *
         diff_in_est * trans(diff_in_est) - (mat1 + trans(mat1)) / scal2;
}


// [[Rcpp::export]]
double BacktrackLineSearch(double kAlpha, double kBeta,
                           NumericVector current_est,
                           NumericVector current_direction,
                           NumericVector current_derv,
                           NumericMatrix u,
                           NumericVector ubar,
                           arma::vec treat, double theta) {

  // A simple backtracking line search
  // This formulation is taken from Algorithm 9.2 of
  // Boyd, Stephen, and Lieven Vandenberghe. Convex optimization. 2004.
  // Args:
  //   kAlpha, kBeta: Constant factors for the backtracking algorithm.
  //   current_est: The current estimate or parameter
  //   current_direction: The descent direction vector for our current
  //                      iteration/step.
  //   current_derv: The derivative of our objective function at the
  //                 current estimate (current_est).
  //   u, ubar, treat, theta: See previous functions above.
  //
  // Returns:
  //   A scalar step size for the BFGS algorithm.

  int N = u.nrow();  // N is the number of subjects
  arma::rowvec cur_derv(current_derv.begin(), current_derv.size(), false);
  arma::colvec cur_direction(current_direction.begin(), current_direction.size(), false);

  arma::mat u_mat(u.begin(), N , u.ncol(), false);
  arma::rowvec u_bar(ubar.begin(), ubar.size(), false);

  // Begin with an initial step size of 1
  double step = 1.0;

  // Initialize some quantities which do not change
  double cur_objective = ObjectiveFunction(current_est, u, ubar, treat, theta);
  // Evaluate the crossprod  of derivative and direction
  double dervf_trans_dervx =  as_scalar(cur_derv * cur_direction);

  // Loop until the objective function of the new
  // proposed step is sufficiently small
  while(ObjectiveFunction(current_est + step * current_direction,
                          u, ubar, treat, theta) >
        cur_objective + kAlpha * step * dervf_trans_dervx) {
    // Reduce step size by factor kBeta.
    step = kBeta * step;
  }

  return step;  // Return step size.
}

// [[Rcpp::export]]
List BFGSAlgorithm(NumericVector initial,
                    NumericMatrix u, NumericVector ubar,
                    arma::vec treat, double theta,
                    double kAlpha, double kBeta,
                    int max_iter, double tol) {
  // Args:
  //   initial: Initial estimate/starting point for BFGS.
  //   u, ubar, treat, theta: See previous functions above.
  //   kAlpha, kBeta: Constant factors for the backtracking algorithm.
  //   max_iter: Maximum number of iterations for algorithm.
  //   tol: Tolerance used in fitting algorithm.
  //
  // Returns:
  //   A List with the following objects
  //     result: The resulting estimate.
  //     weights: The weights for the covariate balancing
  //              methods. See package Vignette.
  //     converge: A boolean indicating convergence.

  //Initialize some variables.
  int K = initial.size();  // Dimension of parameter vector.
  int N = u.nrow();  // Number of subjects.
  arma::mat u_mat(u.begin(), N, u.ncol(), false);

  //Initialize variables which will be used in the main loop.
  arma::vec current_est(initial.begin(), initial.size(), false);
  arma::vec new_est;

  // Old and  new derivative of the objective.
  arma::vec current_derv, new_derv;
  // Old and new descent direction vector.
  arma::vec current_direction, new_direction;

  // Old and new 'estimated' inverse Hessian.
  arma::mat current_inv_hessian = eye<mat>(K, K);
  arma::mat new_inv_hessian;

  //Step size for BFGS and value of objective function.
  double step_size, objective_value;

  // Quantities used for updating inverse Hessian.
  arma::vec diff_in_derv, diff_in_est;

  for (int i = 1; i <= max_iter; i++) {
    // Calculate the current derivative of the objective function and direction.
    current_derv = trans(ObjectiveFirstDerivative(as<NumericVector>(wrap(current_est)),
                                                  u, ubar, treat, theta));
    current_direction = -current_inv_hessian * current_derv;

    // Evaluate the objective function as the current point,
    // we use this later to check if the objective is unbounded.
    objective_value = ObjectiveFunction(as<NumericVector>(wrap(current_est)),
                                        u, ubar, treat, theta);

    //Calculate the Stepsize using the backtracking line search algorithm
    step_size = BacktrackLineSearch(kAlpha, kBeta,
                                    as<NumericVector>(wrap(current_est)),
                                    as<NumericVector>(wrap(current_direction)),
                                    as<NumericVector>(wrap(current_derv)),
                                    u, ubar, treat, theta);

    //Calulate the new estimate
    new_est = current_est + step_size * current_direction;

    //Update the estimate for the inverse hessian using the BFGS update
    new_derv = trans(ObjectiveFirstDerivative(as<NumericVector>(wrap(new_est)),
                                              u, ubar, treat, theta));
    diff_in_derv = new_derv - current_derv;
    diff_in_est = new_est - current_est;
    new_inv_hessian = UpdateInverseHessian(current_inv_hessian,
                                           diff_in_est, diff_in_derv);

    //For some cases, say theta = -1, the objective function can be unbounded
    //in which case the alorithm will never converge.
    //In this case we throw a warning suggesting using a different theta value.
    if (objective_value < -1e+30) {
      // Call the warning function from R.
      Function warning("warning");
      warning("Objective function is unbounded, try different theta values.");
      // Calculate vector of weights.
      arma::vec weights = zeros<vec>(N);
      arma::uvec indx = find(treat == 1);
      arma::vec lambda(new_est.begin(), new_est.size(), false);
      weights.elem(indx) = (CRFamilyDerivative( u_mat * lambda, theta)).elem(indx) / N;
      return List::create(Named("results") = new_est,
                          Named("weights") = weights,
                          Named("converge") = false);
    }

    //Check convergence by first derivative condition.
    if (sum(square(new_derv)) < tol ) {
      //Calculate the weights.
      arma::vec weights = zeros<vec>(N);
      arma::uvec indx = find(treat == 1);
      arma::vec lambda(new_est.begin(), new_est.size(), false);
      weights.elem(indx) = CRFamilyDerivative(u_mat * lambda,
                                              theta).elem(indx) / N;
      return List::create(Named("results") = new_est,
                          Named("weights") = weights,
                          Named("converge") = true);
    } else {
      //If stopping condition is not met update the current estimate and
      //current Inverse Hessian
      current_est = new_est;
      current_inv_hessian = new_inv_hessian;
    }
  }

  //If algotihm runs for maximum number of iterations
  //return the last result and state convergence status is false
  arma::vec weights = zeros<vec>(N);
  arma::uvec indx = find(treat == 1);
  arma::vec lambda(new_est.begin(), new_est.size(), false);
  weights.elem(indx) = CRFamilyDerivative(u_mat * lambda, theta).elem(indx) / N;

  return List::create(Named("results") = new_est,
                      Named("weights") = weights,
                      Named("converge") = false);
}
