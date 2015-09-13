
GetEESimple <- function(u.mat.lambda.treat, u.mat.lambda.placebo,
                        tau.one, tau.zero, Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of simple binary treatments.
  # Args:
  #   u.mat.lambda.treat: The vector u.mat * lambda.treat.
  #   u.mat.lambda.placebo: The vector u.mat * lambda.placebo.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambda.treat and lambda.placebo.

  # gk1 is the EE for lambda_treat (lambda_p in vignette).
  temp1 <- treat * CRFamilyDerivative(u.mat.lambda.treat, theta)
  gk1 <- apply(u.mat, 2, "*", temp1) - u.mat

  # The EE for lambda_placebo (lambda_q in vignette).
  temp2 <- (1 - treat) * CRFamilyDerivative(u.mat.lambda.placebo, theta)
  gk2 <- apply(u.mat, 2, "*", temp2) - u.mat

  # The EE for tau.one, the estimate for E[Y(1)]
  gk3 <- (treat * Y) * CRFamilyDerivative(u.mat.lambda.treat, theta) - tau.one
  # The EE for tau.zero, the estimate for E[Y(0)]
  gk4 <- ((1 - treat) * Y) *
         CRFamilyDerivative(u.mat.lambda.placebo, theta) -
         tau.zero

  # Return the full set of estimatin equations
  cbind(gk1, gk2, gk3, gk4)
}

GetCovSimple <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameters taus.
  # I.e. cov of E[Y(1)] and E[Y(0)]
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstSimple.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The large covariance matrix for the taus,
  #   I.e. cov matrix for tau.one and tau.zero.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  N <- length(Y)  # Sample size.
  n1 <- sum(treat)  # No. of subjects in treatment arm.
  K<- ncol(u.mat)  # No. of covariates.

  # Obtain the lambda and u * lambda objects
  lambda.treat <- (fitted.point.est$lambda.treat)
  lambda.placebo <- (fitted.point.est$lambda.placebo)
  u.mat.lambda.treat <- u.mat %*% lambda.treat
  u.mat.lambda.placebo <- u.mat %*% lambda.placebo

  # Our point estimates tau's
  tau.one <- fitted.point.est$tau.one
  tau.zero <- fitted.point.est$tau.zero

  # Obtain the estimating equations
  gk <- GetEESimple(u.mat.lambda.treat, u.mat.lambda.placebo,
                    tau.one, tau.zero, Y, treat, u.mat, theta)

  # Calculate the MEAT for the Huber-White sandwich estimate
  meat <- crossprod(gk)/N

  # To calculate the bread of the sanwich estimator
  # we follow the notation of the package Vignette.
  # First we create the 2K*2K matrix A
  temp1 <- treat * CRFamilySecondDerivative(u.mat.lambda.treat, theta)
  temp2 <- (1 - treat) * CRFamilySecondDerivative(u.mat.lambda.placebo, theta)

  # We define the matrix A
  tempA1 <- apply(u.mat, 2, "*", temp1)
  A1 <- crossprod(tempA1, u.mat) / N
  tempA2 <- apply(u.mat, 2, "*", temp2)
  A2 <- crossprod(tempA2, u.mat) / N
  A <- bdiag(A1, A2)

  # Now for the 2*2K matrix C
  C1 <- crossprod(u.mat, temp1 * Y)
  C2 <- crossprod(u.mat, temp2 * Y)
  C <- t(bdiag(C1, C2)) / N

  # Calculate the bread of the sandwich estimator
  # using block matrix inverse formula
  Ainv <- bdiag(solve(A1), solve(A2))
  tempMat <- Matrix(0, ncol = 2, nrow = 2 * K)
  bread <- cbind(rbind(Ainv, C %*% Ainv), rbind(tempMat, -diag(rep(1, 2))))

  # Obtain the large covariance matrix
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N

  # Return the submatrix of cov for tau1 and tau0.
  large.cov.mat[-(1:(2 * K)), -(1:(2 * K))]
}

################################################

GetEEATT <- function(u.mat.lambda.placebo, tau.one,
                     tau.zero, Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of estimating the treatment effect on the treated.
  # Args:
  #   u.mat.lambda.placebo: The vector u.mat * lambda.placebo.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambda.placebo.

  n1 <- sum(treat)  # No. in treatment arm.
  N <- length(Y)  # Sample size.
  delta <- n1/N  # Parameter delta.

  #The EE for lambda.placebo
  temp1 <- (1 - treat) * CRFamilyDerivative(u.mat.lambda.placebo, theta)

  gk1 <- apply(u.mat, 2, "*", temp1) - apply(u.mat, 2, "*", treat) / delta
  # THe EE for delta = E[T]
  gk2 <- treat - delta
  # The EE for expected value of tau.one = E[Y(1)|T=1]
  gk3 <- treat * (Y - tau.one) / delta
  # THe EE for expected value of tau.zero = E[Y(0)|T=1]
  gk4 <- temp1 * Y - tau.zero

  cbind(gk1, gk2, gk3, gk4)
}


GetCovATT <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameter vector.
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstATT.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The covariance matrix for our parameters
  #   tau.one and tau.zero.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  # Initialize some variables.
  N <- length(Y)
  n1 <- sum(treat)
  # Recall: delta = E[T], probability of treatment assignment.
  delta <- n1/N
  K<- ncol(u.mat)

  lambda.placebo <- (fitted.point.est$lambda.placebo)
  u.mat.lambda.placebo <- u.mat %*% lambda.placebo

  tau.one <- fitted.point.est$tau.one
  tau.zero <- fitted.point.est$tau.zero

  # Obtain EE and meat for sandwich estimator
  gk <- GetEEATT(u.mat.lambda.placebo = u.mat.lambda.placebo,
                 tau.one = tau.one,
                 tau.zero = tau.zero, Y = Y,
                 treat = treat, u.mat = u.mat,
                 theta = theta)
  meat <- crossprod(gk) / N

  # Obatain the (K+1)*(K+1) matrix A
  temp1 <- (1 - treat) * CRFamilySecondDerivative(u.mat.lambda.placebo, theta)
  tempA <- apply(u.mat, 2, "*", temp1)
  A <- cbind(crossprod(tempA, u.mat), crossprod(u.mat, treat)/(delta ^ 2)) / N
  A <- rbind(A, c(rep(0, K), -1))

  # Evaluate matrix C of dimension 2*(K+1)
  C <- Matrix(0, ncol = K, nrow = 2)
  C[2, ] <- crossprod(u.mat, temp1 * Y)
  C <- cbind(C, c(0, 0)) / N

  # Evaluate the bread matrix and find covariance est
  A.inv <- solve(A)
  temp.mat <- Matrix(0, ncol = 2, nrow = K + 1)
  bread <- cbind(rbind(A.inv, C %*% A.inv), rbind(temp.mat, -diag(rep(1, 2))))

  # Obtain large covariance matrix.
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N
  # Return the submatrix for covariance of tau0 and tau1.
  large.cov.mat[-(1:(K+1)), -(1:(K+1))]
}

################################################

GetEEMultiple <- function(u.mat.lambdas, lambdas, taus,
                          Y, treat, u.mat, theta) {
  # This function calculates the estimating equations for
  # the case of multiple treatment arms.
  # Args:
  #   u.mat.lambdas: The matrix u.mat * lam.mat.
  #   ...: Other parameters to be passed to the function.
  # Returns:
  #   A matrix with N rows, each row is the EE for the
  #   corresponding subject evaluated at the given values of
  #   lambdas.

  # Get the number of treatment arms
  J <- nrow(lambdas)
  # gk1 is a list of EE for lambda_0, lambda_1, ...
  gk1.list <- vector("list", J)
  # gk2 is the list of EE for tau_0, tau_1, ...
  gk2.list <- vector("list", J)

  # For each treatment arm we obtain the EE for
  # lambda_j and tau_j = E[Y(j)]
  for (j in 1:J) {
    u.mat.lambda.j <- u.mat.lambdas[, j]
    temp1 <- 1 * (treat == j - 1) * CRFamilyDerivative(u.mat.lambda.j, theta)
    gk1.list[[j]] <- apply(u.mat, 2, "*", temp1) - u.mat
    gk2.list[[j]] <- 1 * (treat == j - 1) *
                    Y * CRFamilyDerivative(u.mat.lambda.j, theta) -
                    taus[j]
  }

  # Combine the EE for each treatment arm.
  cbind(do.call(cbind, gk1.list), do.call(cbind, gk2.list))
}


GetCovMultiple <- function(fitted.point.est, ...) {
  # This function calculates the estimate of the
  # covariance matrix for the parameter vector with
  # multiple treatments.
  # The function returns the J * J covariance matrix
  # for the parameters tau_0, tau_1, ..., tau_{J-1}
  # Args:
  #   fitted.point.est: The output list of the function
  #                     GetPointEstMultiple.
  #   ...: Other parameters to be passed.
  # Returns:
  #   The covariance matrix for our parameters
  #   tau0, tau1, ..., tau_{J-1}.

  #Obtain extra argments
  args <- list(...)
  Y <- args$Y
  treat <- args$treat
  u.mat <- args$u.mat
  theta <- args$theta

  # Intialize some variables.
  N <- length(Y)
  K <- ncol(fitted.point.est$lam.mat)
  J <- length(unique(treat))

  # Obtain the estimated lambdas and taus.
  lambdas <- fitted.point.est$lam.mat
  taus <- fitted.point.est$tau.treatment.j
  # Obtain the u.mat * lambdas matrix.
  u.mat.lambdas <- tcrossprod(u.mat, lambdas)

  # The meat matrix as usual.
  gk <- GetEEMultiple(u.mat.lambdas = u.mat.lambdas, lambdas = lambdas,
                      taus = taus, Y = Y, treat = treat,
                      u.mat = u.mat,
                      theta = theta)
  meat <- crossprod(gk) / N

  # Alist is the list of length J.
  # Each element is the K * K matrix, a subset of
  # the large A JK * JK matrix.
  A.list <- vector("list", J)
  # This list stores the inverse of each element of A.list.
  A.inv.list <- vector("list", J)

  # And also the C matrix of size J * JK.
  C.list <- vector("list", J)
  for (j in 1:J) {
    u.mat.lambda.j <- u.mat.lambdas[, j]
    temp1 <- 1 * (treat == j - 1) *
             CRFamilySecondDerivative(u.mat.lambda.j, theta)

    #print(dim(u.mat.lambdas))
    # We define the matrix A.
    tempA <- apply(u.mat, 2, "*", temp1)


    A.list[[j]] <- crossprod(tempA, u.mat) / N
    A.inv.list[[j]] <- solve(A.list[[j]])
    C.list[[j]] <- t(crossprod(u.mat, temp1 * Y) / N)
  }

  A.inv <- bdiag(A.inv.list)
  C <- bdiag(C.list)
  temp.mat <- Matrix(0, ncol = J, nrow = J * K)

  # Obtain the bread matrix.
  bread <- cbind(rbind(A.inv, C %*% A.inv),
                 rbind(temp.mat, -diag(rep(1, J))))

  # Obtain full covariance matrix.
  large.cov.mat <- (bread %*% tcrossprod(meat, bread)) / N
  # Return the sub matrix for covariance of Taus.
  large.cov.mat[-(1:(J * K)), -(1:(J * K))]
}
