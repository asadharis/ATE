
GetPointEstSimple <- function(initial.one, initial.two, ...) {
  # Function to get point estimates for the simple case of
  # binary treatment and NOT the case of treatement effect on
  # the treated.
  #
  # Args:
  #   initial.one: The initial vector for first BFGS algorithm (treatment arm).
  #   initial.two: The initial vector for second BFGS algorithm (placebo arm).
  #   ...: Other arguments to be passed to the function from parent.
  #
  # Returns:
  #   A list of estimated lambda values where lambda is the parameter
  #   vector which we optiize over. It also contains the weights
  #   obtained from the estimated lambda values used for covariate
  #   balancing estimates. Finally, the list also contains the
  #   point estimates and convergence indicator.

  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  u.mat <- cbind(1, X)  # Obtain u matrix of covariates.
  u.bar <- colMeans(u.mat)  # Get column means.

  # Perform some simple checks.
  if (length(initial.one) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (length(initial.two) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("Running BFGS algorithm for estimating Weights p: \n")
  }

  # Implement the BFGS algorithm to get the optimal lambda
  # and correspoding weights.
  # This corresponds to \lambda_p and \hat{p}_K in vignette.
  treatment.hat <- BFGSAlgorithm(initial.one, u.mat, u.bar, treat,
                                 theta, backtrack.alpha,
                                 backtrack.beta, max.iter,
                                 tol)

  if (verbose) {
    cat("\nRunning BFGS algorithm for estimating Weights q: \n")
  }

  # Implement the BFGS algorithm to get the optimal lambda_q
  # and correspoding weights \hat{q}_K.
  placebo.hat <- BFGSAlgorithm(initial.two, u.mat, u.bar, 1 - treat,
                               theta, backtrack.alpha,
                               backtrack.beta,
                               max.iter, tol)

  # Obtain estimates for tau1 = E[Y(1)], tau2 = E[Y(0)]
  # and tau = E[Y(1)] - E[Y(0)].
  tau.one  <- sum((treatment.hat$weights * Y)[treat == 1])
  tau.zero <- sum((placebo.hat$weights * Y)[treat == 0])
  tau <- tau.one - tau.zero

  # Obtain estimated weights and lambda's for
  # each treatment arm.
  weights.treat <- treatment.hat$weights
  weights.placebo <- placebo.hat$weights

  lambda.treat <- treatment.hat$results
  lambda.placebo <- placebo.hat$results

  # Throw a warning if algorithm did not converge.
  converge = TRUE
  if (!treatment.hat$converge || !placebo.hat$converge) {
    warning("BFGS Algorithm did not converge for atleast one objective function.")
    converge <- FALSE
  }

  # Return list of objects
  return(list(lambda.treat    = lambda.treat,
              lambda.placebo  = lambda.placebo,
              weights.treat   = weights.treat,
              weights.placebo = weights.placebo,
              tau.one = tau.one, tau.zero = tau.zero,
              tau = tau, converge = converge))
}


GetPointEstATT <- function(initial, ...) {
  # Function to get point estimates average treatment effect
  # on the treated. This case also has binary treatment.
  #
  # Args:
  #   initial: The initial vector for the BFGS algorithm.
  #   ...: Other arguments to be passed to the function.
  #
  # Returns:
  #   List of objects similar to the previous function. However, this
  #   function does not return lambda or weights for the treatment arm.
  #   This is because the we do not need a covariate balancing technique
  #   for the treatment arm in this case.

  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  u.mat <- cbind(1, X)  # Design matrix u.

  # Main difference here is the definition of u.bar.
  # u.bar is vector of column means ONLY for those
  # who recieved the treatment.
  u.bar <- colMeans(u.mat[treat == 1, ])

  if (length(initial) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("\nRunning BFGS algorithm Raphson for estimating Weights q: \n")
  }

  # Run the BFGS algorithm for obtaining the parameter and weights
  # used for covariate balancing.
  placebo.hat <- BFGSAlgorithm(initial, u.mat, u.bar,
                                1 - treat, theta,
                                backtrack.alpha, backtrack.beta,
                                max.iter, tol)

  # Note that treatment effect on treated is simple to estimate
  tau.one <- mean(Y[treat == 1])
  # The calibration estimator for E[Y(0)|T = 1]
  tau.zero <- sum((placebo.hat$weights * Y)[treat == 0])
  tau <- tau.one - tau.zero

  # Weights and lambda vectors
  weights.placebo <- placebo.hat$weights
  lambda.placebo <- placebo.hat$results

  # Warning message if the algorithm did not converge
  converge = TRUE
  if (!placebo.hat$converge) {
    warning("\nBFGS algorithm did not converge for the objective function.")
    converge <- FALSE
  }

  # Return list
  return(list(lambda.placebo = lambda.placebo,
              weights.placebo = weights.placebo,
              tau.one = tau.one, tau.zero = tau.zero,
              tau = tau, converge = converge))
}


GetPointEstMultiple <- function(initial.mat, ...) {
  # Function to get point estimates for treatment effect
  # when we have multiple treatment arms.
  #
  # Args:
  #   initial.mat: A matrix of initial values for the different
  #                BFGS algorithm. Each row is an inital vector for
  #                an algorithm. The total number of rows is J: number of
  #                different treatment arms.
  #   ...: Other arguments to be passed to the function.
  #
  # Returns:
  #   A List with the following objects
  #     lam.mat: A matrix of estimated lambda values. This has the same
  #              dimensions as initial.mat.
  #     weights.mat: A matrix estimated weights. The j-th row
  #                  corresponds to the weights for treatment j.
  #     tau.treatment.j: A vector of estimates of [EY(0), EY(1),..., EY(J-1)].
  #     converge: A boolean indicator of convergence status.


  #Obtain extra arguments
  args<- list(...)
  Y <- args$Y
  X <- args$X
  treat <- args$treat
  theta <- args$theta
  max.iter <- args$max.iter
  tol <- args$tol
  backtrack.alpha <- args$backtrack.alpha
  backtrack.beta <- args$backtrack.beta
  verbose <- args$verbose

  N <- length(Y)  # Sample size.
  J <- length(unique(treat))  # Number of treatment arms.
  u.mat <- cbind(1, X)  # Obtain design matrix u.
  u.bar <- colMeans(u.mat)  # Obtain column means of design matrix.

  # A simple verification of dimensions.
  if (ncol(initial.mat) != length(u.mat[1, ])) {
    stop("Incorrect length of initial vector")
  }
  if (verbose) {
    cat("\nRunning BFGS for estimating Weights: \n")
  }

  # Initialize the matrices and vector which will be returned
  # by this function.
  lam.mat <- matrix(0, ncol = ncol(initial.mat), nrow = J)
  weights.mat <- matrix(0, ncol = N, nrow = J)
  tau.treatment.j <- numeric(J)

  # Loop through the different treatment arms.
  for (j in 0:(J - 1)) {
    # Inidicator of treatment arm.
    temp.treat <- 1 * (treat == j)

    # Implement BFGS algorithm
    treatment.j.hat <- BFGSAlgorithm(initial.mat[j + 1, ],
                                     u.mat, u.bar, temp.treat,
                                     theta,
                                     backtrack.alpha,
                                     backtrack.beta, max.iter, tol)

    # Find estimate for E[Y(j)]
    tau.j.hat <- sum((treatment.j.hat$weights * Y)[temp.treat == 1])

    tau.treatment.j[j + 1] <- tau.j.hat
    lam.mat[j + 1, ] <- treatment.j.hat$results
    weights.mat[j + 1, ] <- treatment.j.hat$weights

    # Warning for non-convergence
    converge = TRUE
    if (!treatment.j.hat$converge) {
      warning(paste("BFGS algorithm did not converge for treatment arm",
                    j))
      converge <- FALSE
    }
  }
  # Return result.

  return(list(lam.mat = lam.mat,
              weights.mat = weights.mat,
              tau.treatment.j = tau.treatment.j,
              converge = converge))
}
