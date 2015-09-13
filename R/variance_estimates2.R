
get.gk.simple <- function(u, Y, Ti, lam, beta, tau1, tau0, theta, ulam, ubeta) {
  # gk1 is the EE for lambda_p
  temp1 <- Ti * cpp_dcr_rho(ulam, theta)
  gk1 <- apply(u, 2, "*", temp1) - u

  # The EE for lambda_q
  temp2 <- (1 - Ti) * cpp_dcr_rho(ubeta, theta)
  gk2 <- apply(u, 2, "*", temp2) - u

  # The EE for tau1, the estimate for E[Y(1)]
  gk3 <- (Ti * Y) * cpp_dcr_rho(ulam, theta) - tau1
  # The EE for tau0, the estimate for E[Y(0)]
  gk4 <- ((1 - Ti) * Y) * cpp_dcr_rho(ubeta, theta) - tau0

  cbind(gk1, gk2, gk3, gk4)
}


# This function obtains the big covariance matrix for the parameter theta of Section 2.4
get.cov.simple <- function(u, Y, Ti, obj, theta) {
  # Initialize some objects
  N <- length(Y)
  n1 <- sum(Ti)
  lam <- (obj$lam.p)
  beta <- (obj$lam.q)
  ulam <- u %*% lam
  ubeta <- u %*% beta

  tau1 <- obj$Y1
  tau0 <- obj$Y0
  K <- length(lam)

  # Obtain the estimating equations
  gk <- get.gk.simple(u, Y, Ti, lam, beta, tau1, tau0, theta, ulam, ubeta)
  # Calculate the MEAT for the Huber-White sandwich estimate
  meat <- crossprod(gk)/N

  # Again, we follow the notation of the package Vignette First we create the 2K*2K matrix A
  temp1 <- Ti * cpp_ddcr_rho(ulam, theta)
  temp2 <- (1 - Ti) * cpp_ddcr_rho(ubeta, theta)
  # We define the matrix A
  tempA1 <- apply(u, 2, "*", temp1)
  A1 <- crossprod(tempA1, u)/N
  tempA2 <- apply(u, 2, "*", temp2)
  A2 <- crossprod(tempA2, u)/N
  A <- bdiag(A1, A2)

  # Now for the 2*2K matrix C
  C1 <- crossprod(u, temp1 * Y)
  C2 <- crossprod(u, temp2 * Y)
  C <- t(bdiag(C1, C2))/N

  # Calculate the bread of the sandwich estimator using block matrix inverse formula
  Ainv <- bdiag(solve(A1), solve(A2))
  tempMat <- Matrix(0, ncol = 2, nrow = 2 * K)
  bread <- cbind(rbind(Ainv, C %*% Ainv), rbind(tempMat, -diag(rep(1, 2))))

  # Finally return the covariance matrix
  (bread %*% tcrossprod(meat, bread))/N
}


###############################################################################

# The EE equations for estimate the average treatment effect on the treated. ATT
get.gk.ATT <- function(u, Y, Ti, beta, tau1, tau0, theta, ubeta) {

  n1 <- sum(Ti)
  N <- length(Y)
  frac <- n1/N

  # The notation is slightly mixed here compared to package Vignette
  temp1 <- (1 - Ti) * cpp_dcr_rho(ubeta, theta)
  gk1 <- apply(u, 2, "*", temp1) - apply(u, 2, "*", Ti)/frac
  # THe EE for expected value of E[Y(1)|T=1]
  gk2 <- Ti * (Y - tau1)/frac
  # THe EE for expected value of E[Y(0)|T=1]
  gk3 <- temp1 * Y - tau0
  # THe EE for expected value E[T]
  gk4 <- Ti - frac

  # This fixes the ordering to agree with package vignette
  cbind(gk1, gk4, gk2, gk3)
}


# This function obtains the big covariance matrix for all coefficients for ATT
get.cov.ATT <- function(u, Y, Ti, obj, theta) {
  # We largely follow the same format as above
  N <- length(Y)
  n1 <- sum(Ti)
  frac <- n1/N

  beta <- (obj$lam.q)
  ubeta <- u %*% beta

  tau1 <- obj$Y1
  tau0 <- obj$Y0
  K <- length(beta)

  # Obtain gk and meat for sandwich estimator
  gk <- get.gk.ATT(u, Y, Ti, beta, tau1, tau0, theta, ubeta)
  meat <- crossprod(gk)/N

  # Obatain the (K+1)*(K+1) matrix A
  temp1 <- (1 - Ti) * cpp_ddcr_rho(ubeta, theta)
  tempA <- apply(u, 2, "*", temp1)
  A <- cbind(crossprod(tempA, u), crossprod(u, Ti)/(frac^2))/N
  A <- rbind(A, c(rep(0, K), -1))

  # Evaluate matrix C of dimension 2*(K+1)
  C <- Matrix(0, ncol = K, nrow = 2)
  C[2, ] <- crossprod(u, temp1 * Y)
  C <- cbind(C, c(0, 0))/N

  # Evaluate the bread matrix and find covariance est
  Ainv <- solve(A)
  tempMat <- Matrix(0, ncol = 2, nrow = K + 1)
  bread <- cbind(rbind(Ainv, C %*% Ainv), rbind(tempMat, -diag(rep(1, 2))))

  (bread %*% tcrossprod(meat, bread))/N
}


###############################################################################

# The EE equations for multiple treatment arms
get.gk.MT <- function(u, Y, Ti, lams, taus, theta, ulams) {
  J <- nrow(lams)
  gk1list <- vector("list", J)
  gk2list <- vector("list", J)

  # For each treatment arm we obtain the EE for lamba_j and tau_j = E[Y(j)]
  for (j in 1:J) {
    ulam <- ulams[, j]
    temp1 <- 1 * (Ti == j - 1) * cpp_dcr_rho(ulam, theta)
    gk1list[[j]] <- apply(u, 2, "*", temp1) - u
    gk2list[[j]] <- 1 * (Ti == j - 1) * Y * cpp_dcr_rho(ulam, theta) - taus[j]
  }
  cbind(do.call(cbind, gk1list), do.call(cbind, gk2list))
}

# This function obtains the BIG covariance matrix for all coefficients for MT
get.cov.MT <- function(u, Y, Ti, obj, theta) {
  N <- length(Y)
  lams <- obj$lam.mat
  K <- ncol(lams)
  J <- length(unique(Ti))
  taus <- obj$Yj.hat
  ulams <- tcrossprod(u, lams)

  # The meat matrix as usual
  gk <- get.gk.MT(u, Y, Ti, lams, taus, theta, ulams)
  meat <- crossprod(gk)/N

  # Evaluate the A matrix of size JK*JK
  Alist <- vector("list", J)
  AinvList <- vector("list", J)

  # And also the C matrix of size J*JK
  Clist <- vector("list", J)
  for (j in 1:J) {
    ulam <- ulams[, j]
    temp1 <- 1 * (Ti == j - 1) * cpp_ddcr_rho(ulam, theta)
    # We define the matrix A
    tempA <- apply(u, 2, "*", temp1)
    Alist[[j]] <- crossprod(tempA, u)/N
    AinvList[[j]] <- solve(Alist[[j]])
    Clist[[j]] <- t(crossprod(u, temp1 * Y)/N)

  }

  Ainv <- bdiag(AinvList)
  C <- bdiag(Clist)
  tempMat <- Matrix(0, ncol = J, nrow = J * K)

  # Obtain the bread matrix
  bread <- cbind(rbind(Ainv, C %*% Ainv), rbind(tempMat, -diag(rep(1, J))))
  (bread %*% tcrossprod(meat, bread))/N
}

###############################################################################

# This function takes in an object given by the functions get.est.*
estimate_variance <- function(object) {
  K <- object$K
  u <- cbind(1, object$X)
  if (object$gp == "simple") {
    # get the full covariance matrix
    cov.mat <- get.cov.simple(u, object$Y, object$Ti, object, object$theta)
    # Obtain the lower right corner 2*2 matrix for tau1 and tau0
    covEY <- cov.mat[-(1:(2 * K)), -(1:(2 * K))]
    fin <- covEY

  } else if (object$gp == "ATT") {
    # Again, we get the full cov matrix
    cov.mat <- get.cov.ATT(u, object$Y, object$Ti, object, object$theta)

    # Again the lower right corner is for our relevant parameters
    covEY <- cov.mat[-(1:K), -(1:K)]
    covEY <- covEY[-1, -1]
    fin <- covEY
  } else {
    # Full covariance matrix
    cov.mat <- get.cov.MT(u, object$Y, object$Ti, object, object$theta)
    # The lower right J*J submatrix is for our parameters
    covEY <- cov.mat[-(1:(object$J * K)), -(1:(object$J * K))]
    fin <- covEY
  }
  return(fin)
}
