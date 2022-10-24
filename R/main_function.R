
ATE <- function(Y, treat, X, theta = 1, ATT = FALSE,
                verbose = FALSE, max.iter = 100,
                tol = 1e-10, initial.values = NULL,
                backtrack.alpha = 0.3, backtrack.beta = 0.5) {

  # This is the main function available to the user.
  # This creates an ATE object. The object contains point
  # estimates and variance estimates. It also has
  # generic S3 methods such as print, plot and summary.
  # Args:
  #   Y: Response vector.
  #   X: Design/Covariate matrix.
  #   treat: Vector of treatment indicator. Must be of
  #          the form (0, 1, ...).
  #   ATT: Indicate if treatment effect on the treated is
  #        to be estimated.
  #   verbose: Indicate if extra statments should be printed
  #            to show progress of the function while it runs.
  #   max.iter: Maximum no. of iterations for BFGS algorithms.
  #   tol: Tolerance of the algorithm for stopping conditions.
  #   initial.values: Matrix of initial values for BFGS algorithm.
  #   backtrack.alpha, backtrack.beta: Parameters for backtrack line
  #                                    search algorithm.
  # Returns:
  #   An object of class 'ATE'. This contains point estimates and
  #   and variance covariance matrix of our estimates. Along with other
  #   information about our object such as data used for estimation.

  # J is the number treatment arms
  J <- length(unique(treat))

  # In some cases if we have a data.frame we need
  # to convert it into a numeric matrix
  X_original <- X
  if (is.data.frame(X)) {
    if (any(vapply(X, function(x) !is.numeric(x), logical(1L)))) {
      for (i in names(X)) {
        if (is.character(X[[i]])) {
          X[[i]] <- factor(X[[i]])
          X_original[[i]] <- X[[i]]
        }
        if (is.factor(X[[i]])) {
          new <- model.matrix(reformulate(i), data = X)[,-1, drop = FALSE]
          if (i == names(X)[1L]) {
            X <- cbind(as.data.frame(new), X[,-1,drop = FALSE])
          }
          else if (i == names(X)[ncol(X)]) {
            X <- cbind(X[,-ncol(X), drop = FALSE], as.data.frame(new))
          }
          else {
            X <- cbind(X[,1:(which(names(X) == i) - 1)],
                       as.data.frame(new),
                       X[,(which(names(X) == i) + 1):ncol(X), drop = FALSE])
          }
        }
        else {
          X[[i]] <- as.numeric(X[[i]])
        }
      }
    }
    X <- as.matrix(X)
  }
  # Take care of the case of a single covariate
  else if (is.numeric(X) && length(dim(X)) != 2L) {
    warning("Data matrix 'X' is a vector, will be treated as n x 1 matrix.",
            immediate. = TRUE)
    X <- matrix(X, ncol = 1, nrow = length(X))
    X_original <- as.data.frame(X)
  }
  else if (is.matrix(X)) {
    X_original <- as.data.frame(X)
  }
  else {
    stop("Data matrix 'X' must be a data.frame or matrix.")
  }
  # The total length of the parameter vector Recall:
  #Using the notation of the Package Vignette we have u_K(X) = cbind(1,X).
  K <- ncol(X) + 1

  # Some simple checks before running the main functions
  if (nrow(X) != length(Y)) {
    stop("Dimensions of covariates and response do not match.")
  }
  if (J == 1) {
    stop("There must be atleast two treatment arms")
  }
  if (!all(unique(treat) %in% 0:(J - 1))) {
    stop("The treatment levels must be labelled 0,1,2,...")
  }
  if (!is.numeric(theta)) {
    stop("'theta' must be a real number.")
  }
  if (J > 2 & ATT) {
    stop("For ATT == TRUE, must have only 2 treatment arms.")
  }

  # If the initial values are not
  # provided we generate simple initial parameters.
  if (is.null(initial.values)) {
    if (ATT) {
      initial.values <- numeric(K)
    } else {
      initial.values <- matrix(0, ncol = K, nrow = J)
    }
  }

  # For ATT there is only one implementation of BFGS.
  if (ATT & !is.vector(initial.values)) {
    stop("For ATT == TRUE, only need one vector of initial values.")
  }
  # In all other cases this has to be a matrix.
  if (!ATT) {
    if (!is.matrix(initial.values)) {
      stop("Initial values must be a matrix")
    }
    if (any(dim(initial.values) != c(J, K))) {
      stop("Matrix of initial values must have dimensions J x K.")
    }
  }

  # Now we specify the category of the problem The simple case with binary treatment.
  gp <- "simple"

  # The case of average treatment effect on the treated.
  if (ATT) {
    gp <- "ATT"
  }
  # The case of Multiple treatment arms.
  if (J > 2) {
    gp <- "MT"
  }

  # Obtain the estimates for whatever the case may be.
  if (gp == "simple") {
    est <- GetPointEstSimple(initial.values[1, ],
                             initial.values[2, ],
                             X = X, Y = Y, treat = treat,
                             theta = theta, max.iter = max.iter,
                             tol = tol,
                             backtrack.alpha = backtrack.alpha,
                             backtrack.beta = backtrack.beta,
                             verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovSimple(est, Y = Y, treat = treat,
                            u.mat = cbind(1, X),
                            theta = theta)

  } else if (gp == "ATT") {
    est <- GetPointEstATT(initial.values, X = X, Y = Y,
                          treat = treat,
                          theta = theta, max.iter = max.iter,
                          tol = tol,
                          backtrack.alpha = backtrack.alpha,
                          backtrack.beta = backtrack.beta,
                          verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovATT(est, Y = Y, treat = treat,
                         u.mat = cbind(1, X),
                         theta = theta)
  } else if (gp == "MT") {
    est <- GetPointEstMultiple(initial.values, X = X, Y = Y,
                               treat = treat,
                               theta = theta, max.iter = max.iter,
                               tol = tol,
                               backtrack.alpha = backtrack.alpha,
                               backtrack.beta = backtrack.beta,
                               verbose = verbose)
    if (verbose) {
      cat("\nEstimating Variance")
    }
    cov.mat <- GetCovMultiple(est, Y = Y, treat = treat,
                              u.mat = cbind(1, X),
                              theta = theta)
  }

  # Begin building the 'ATE' object.
  res <- est
  res$vcov<- cov.mat
  res$X <- X_original
  res$Y <- Y
  res$treat <- treat
  res$theta <- theta
  res$gp <- gp
  res$J <- J
  res$K <- K
  res$call <- match.call()

  # Rename some of the elements of the list and organize the output object
  if (gp == "simple") {
    estimate <- c(res$tau.one, res$tau.zero, res$tau)
    names(estimate) <- c("E[Y(1)]", "E[Y(0)]", "ATE")
    res$estimate <- estimate
    res$tau.one <- NULL
    res$tau.zero <- NULL
    res$tau <- NULL
  } else if (gp == "ATT") {
    estimate <- c(res$tau.one, res$tau.zero, res$tau)
    names(estimate) <- c("E[Y(1)|T=1]", "E[Y(0)|T=1]", "ATT")
    res$estimate <- estimate
    res$tau.one <- NULL
    res$tau.zero <- NULL
    res$tau <- NULL
  } else {
    estimate <- res$tau.treatment.j
    names(estimate) <- paste("E[Y(", 0:(J - 1), ")]", sep = "")
    res$estimate <- estimate
    res$tau.treatment.j <- NULL
  }
  # Define the class ATE
  class(res) <- "ATE"
  return(res)
}

print.ATE <- function(x, ...) {
  # S3 print method for class ATE.
  # The function prints the point estimates
  # Args:
  #   x: Object of class 'ATE'.
  #   ...: Other arguments. Note: these are not used.
  # Returns:
  #   Prints the function call and point estimates.

  if (x$gp == "simple") {
    cat("Call:\n")
    print(x$call)
    cat(paste0("\nThe analysis was completed for a ",
               "simple study design with binary treatment.\n"))
    cat("\nPoint Estimates:\n")
    print(x$estimate)
  } else if (x$gp == "ATT") {
    cat("Call:\n")
    print(x$call)
    cat(paste0("\nThe analysis was completed for a ",
               "binary treatment for estimating treatment ",
               "effect on the treated.\n"))
    cat("\nPoint Estimates:\n")
    print(x$estimate)
  } else {
    cat("Call:\n")
    print(x$call)
    cat(paste0("\nThe analysis was completed for ",
               "a study design with multiple treatment arms.\n"))
    cat("\nPoint Estimates:\n")
    print(x$estimate)
  }
  invisible(x)
}

summary.ATE <- function(object, ...) {
  # S3 summary method for ATE. This function calculates
  # the SE, Z-statistic and confidence intervals
  # and P-values.
  # Args:
  #   object: An object of class 'ATE'.
  #   ...: Other arguments (not used).
  # Returns:
  #   An object of type summary.ATE. This is a list with a matrix for
  #   coefficient estimates. This matrix contains, SE, Z-statistic etc.


  # For binary treatment we calculate the variance of EY(1) - EY(0)
  # using the formula var(a-b) = var(a) + var(b) - 2cov(a,b).
  if (object$gp == "simple" || object$gp == "ATT") {
    var.tau <- object$vcov[1, 1] + object$vcov[2, 2] - 2 * object$vcov[1, 2]
    se <- c(sqrt(diag(object$vcov)), sqrt(var.tau))
  } else {
    # The case of multiple treatments, now we do not have a
    # specific variable like EY(1) - EY(0). So we just obtain
    # the SE for EY(0), EY(1), EY(2),... .
    se <- sqrt(diag(object$vcov))
  }

  # Obtain confidence intervals
  ci.lower <- object$estimate + se * qnorm(0.025)
  ci.upper <- object$estimate + se * qnorm(0.975)

  # Evaluate the Z-statistic and P-value
  z.stat <- object$estimate / se
  p.values <- 2 * pnorm(-abs(z.stat))

  # Create a coefficient matrix which we can print later
  # similar to the out-put of summary for lm/glm.
  coef <- cbind(Estimate = object$estimate,
                'Std. Error' = se,
                '95%.Lower' = ci.lower, '95%.Upper' = ci.upper,
                'z value' = z.stat,
                'p value' = p.values)

  if (object$gp == "simple") {
    res <- list(call = object$call, Estimate = coef,
                vcov = object$vcov, converge = object$converge,
                weights.treat = object$weights.treat,
                weights.placebo = object$weights.placebo)

  } else if (object$gp == "ATT") {
    res <- list(call = object$call, Estimate = coef,
                vcov = object$cov, converge = object$converge,
                weights.treat = object$weights.treat,
                weights.placebo = object$weights.placebo)

  } else {
    res <- list(call = object$call, Estimate = coef,
                vcov = object$cov, converge = object$converge,
                weights = object$weights.mat)
  }
  class(res) <- "summary.ATE"
  return(res)
}

print.summary.ATE <- function(x, ...) {
  # A Print method for 'summary.ATE' which uses the
  # coefficient matrix to give us print statements like the summary
  # function of lm/glm.
  # Args:
  #   x: An object of type 'summary.ATE'.
  #   ...: Other arguments to be passed to the function.
  # Returns:
  #   Prints an output matrix similar to the out put of lm/glm
  #   functions.

  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$Estimate, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

################################################

weighted.ecdf <- function(t, x, weights = NULL) {
  # A simple function to evaluate the emiprical
  # CDF or the  weighted emiprical CDF.
  # Args:
  #   t: A scalar point at which we wish to evaluate the eCDF.
  #   x: The data points used for evaluating the eCDF.
  #   weights: The vector of weights for the weighted CDF.
  # Returns:
  #   Plots the empirical CDF or weighted Empirical CDFs at t.

  if (is.null(weights)) {
    sum(1 * (x <= t)) / length(x)
  } else {
    sum(1 * (x <= t) * weights)
  }
}

plot.ATE <- function(x, ...) {
  # The S3 plot function for ATE.
  # This function plots the empirical CDF for all the
  # covariates for the different
  # tretment groups along with theweighted eCDF.
  # This shows the effect of covariate balancing.
  # Args:
  #   x: An object of class 'ATE'.
  #   ...: Other arguments to be passed.
  # Returns:
  #   Plots for the eCDF and weighted eCDF for each covariate.

  treat <- x$treat

  # Split factor variables into dummies
  if (any(vapply(x$X, function(x) !is.numeric(x), logical(1L)))) {
    x$X <- as.data.frame(
      model.matrix(reformulate(names(x$X)), data = x$X,
                   contrasts.arg = lapply(Filter(is.factor, x$X), contrasts,
                                          contrasts = FALSE))
    )[,-1, drop = FALSE]
  }

  ################## Case 1: Binary treatment###################

  if (x$gp == "simple" || x$gp == "ATT") {
    # Boolean to see if ATT is true.
    ATT <- x$gp == "ATT"

    # Obtain weights for each treatment arm.
    if (!ATT) {
      weights.treat <- x$weights.treat[treat == 1]
    }
    weights.placebo <- x$weights.placebo[treat == 0]

    # Check if covariates are named.
    names <- colnames(x$X)
    # If not, assign generic names.
    if (is.null(names)) {
      names <- paste("X", 1:ncol(x$X), sep = "")
    }

    # Obtain the subsample of subjects for each treatment arm
    x.treat <- x$X[treat == 1,]
    x.placebo <- x$X[treat == 0,]
    for (i in 1:ncol(x.treat)) {
      # Plot the first covariate Allow user to
      # sequentially view subsequent plots
      if (i == 2) {
        par(ask = TRUE)
      }

      # A special plot function for binary covariates
      if (length(unique(x$X[[i]])) == 2) {
        Treatment <- x.treat[[i]]
        Placebo <- x.placebo[[i]]
        # First plot the unweighted case.
        # This function plots the means for each treatment.
        par(mfrow = c(1, 2))
        # Plot dots for the mean for treatment and placebo at
        # position 1 and 2.
        plot(c(0.5, 1, 2, 2.5), c(2, mean(Treatment), mean(Placebo), 2),
             pch = 16, cex = 1.5, ylim = range(x$X[[i]]),
             ylab = "Mean of group",
             xlab = "", col = c("blue", "red"),
             main = paste0("Unweighted\n", names[i]), xaxt = "n")
        axis(side = 1, at = c(1, 2),
             labels = c("Treatment", "Placebo"))
        if (ATT) {
          abline(h = mean(Treatment), lty = 2)
        } else {
          abline(h = mean(x$X[[i]]), lty = 2)
        }
        # Now for the weighed means for treatment and placebo group.
        if (!ATT) {
          new.treat <- sum(weights.treat * Treatment)
        } else {
          new.treat <- mean(Treatment)
        }
        new.placebo <- sum(weights.placebo * Placebo)
        # Again, we plot the means for each treatment group at
        # position 1 and 2.
        plot(c(0.5, 1, 2, 2.5), c(2, new.treat, new.placebo, 2),
             pch = 16, cex = 1.5, ylim = range(x$X[[i]]),
          ylab = "Mean of group", xlab = "",
          col = c("blue", "red"), main = paste0("Weighted\n", names[i]),
          xaxt = "n")
        axis(side = 1, at = c(1, 2), labels = c("Treatment", "Placebo"))
        if (ATT) {
          abline(h = mean(Treatment), lty = 2)
        } else {
          abline(h = mean(x$X[[i]]), lty = 2)
        }

      } else { # Plot for continuous covariates.

        # Obtain the range and initialize a
        # sequence at which we will plot the CDF.
        # my.range <- range(x$X[[i]])
        # my.seq <- seq(my.range[1], my.range[2], length = 100)
        my.seq <- sort(unique(x$X[[i]]))

        ecdf.treat <- sapply(my.seq, weighted.ecdf, x = x.treat[[i]])
        ecdf.placebo <- sapply(my.seq, weighted.ecdf, x = x.placebo[[i]])

        # Plot the unweighted empirical CDF for each case
        par(mfrow = c(1, 2))
        plot(my.seq, ecdf.treat, cex = 0.4, pch = 16,
             type = "s", lty = 1, col = "red",
             xlab = names[i], ylab = "empirical CDF",
             main = "Unweighted empirical CDF")
        lines(my.seq, ecdf.placebo, cex = 0.4,
              type = "s", pch = 16, lty = 2, col = "blue")
        legend("bottomright", c("Treatment", "Placebo"),
               lty = c(1, 2), col = c("red", "blue"))

        # Now we plot the weighted eCDF
        if (!ATT) {
          weighted.ecdf.treat <- sapply(my.seq, weighted.ecdf,
                                        x = x.treat[[i]],
                                        weights = weights.treat)
        } else {
          weighted.ecdf.treat <- ecdf.treat
        }
        weighted.ecdf.placebo <- sapply(my.seq, weighted.ecdf,
                                        x = x.placebo[[i]],
                                        weights = weights.placebo)
        plot(my.seq, weighted.ecdf.treat, cex = 0.4, pch = 16,
             type = "s", lty = 1, col = "red",
             xlab = names[i], ylab = "Empirical CDF",
             main = "Weighted Empirical CDF")
        lines(my.seq, weighted.ecdf.placebo, cex = 0.4,
              type = "s", pch = 16, lty = 2, col = "blue")
      }
    }
    par(ask = FALSE)

  ################## Case 2: Multiple Treatments ###################

  } else {
    weights.mat <- x$weights.mat
    # Check names or assign generic names.
    names <- colnames(x$X)
    if (is.null(names)) {
      names <- paste("X", 1:ncol(x$X), sep = "")
    }

    # We generalize the case of binary treatment to
    # J treatment arms.
    # Much of the syntax is the same.
    J <- x$J
    for (i in 1:ncol(x$X)) {
      if (i == 2) {
        par(ask = TRUE)
      }

      # Again, for binary covariates first.
      if (length(unique(x$X[[i]])) == 2) {

        par(mfrow = c(1, 2))
        # Plot the means for each group at positions 0, 1, ..., J-1.
        plot(c(0:(J - 1)), c(mean(x$X[[i]][treat == 0]), rep(2, J - 1)),
             pch = 16, cex = 1.5, ylim = range(x$X[[i]]),
             ylab = "Mean of group", xlab = "Treatment group",
             col = 1, main = paste0("Unweighted\n", names[i]),
             xaxt = "n",
             xlim = c(-0.5, J - 1 + 0.5))
        axis(side = 1, at = 0:(J - 1), labels = paste(0:(J - 1)))
        abline(h = mean(x$X[[i]]), lty = 2)

        for (j in 1:(J - 1)) {
          points(j, mean(x$X[[i]][treat == j]), pch = 16,
                 cex = 1.5, col = j + 1)
        }

        #Now for the weighted means.
        plot(c(0:(J - 1)),
             c(sum(x$X[[i]][treat == 0] * weights.mat[1, treat == 0]),
             rep(2, J - 1)), pch = 16, cex = 1.5,
             ylim = range(x$X[[i]]), ylab = "Mean of group",
             xlab = "Treatment group", col = 1,
             main = paste0("Weighted\n", names[i]),
             xaxt = "n", xlim = c(-0.5, J - 1 + 0.5))
        axis(side = 1, at = 0:(J - 1), labels = paste(0:(J - 1)))
        abline(h = mean(x$X[[i]]), lty = 2)
        for (j in 1:(x$J - 1)) {
          points(j, sum(x$X[[i]][treat == j] * weights.mat[j + 1, treat == j]),
                 pch = 16, cex = 1.5, col = j + 1)
        }

      } else {  # The case of continuous variables.
        # The Set-up is virtually unchaged.
        # my.range <- range(x$X[[i]])
        # my.seq <- seq(my.range[1], my.range[2], length = 100)
        my.seq <- sort(unique(x$X[[i]]))

        ecdf.placebo <- sapply(my.seq, weighted.ecdf, x = x$X[[i]][x$treat == 0])
        par(mfrow = c(1, 2))
        plot(my.seq, ecdf.placebo, cex = 0.4, pch = 16, type = "s",
             lty = 1, col = 1, xlab = names[i],
             ylab = "Empirical CDF", main = "Unweighted Empirical CDF")

        # A for loop to cycle through all the others
        for (j in 1:(J - 1)) {
          ecdf.treat.j <- sapply(my.seq, weighted.ecdf, x = x$X[[i]][treat == j])
          lines(my.seq, ecdf.treat.j, cex = 0.4, pch = 16,
                type = "s", lty = j + 1, col = j + 1)
        }
        legend("bottomright", paste("gp", 0:(x$J - 1)),
               lty = 1:J, col = 1:J)

        weighted.ecdf.placebo <- sapply(my.seq, weighted.ecdf,
                               x = x$X[[i]][treat == 0],
                               weights = weights.mat[1, treat == 0])
        plot(my.seq, weighted.ecdf.placebo, cex = 0.4, pch = 16,
             type = "s", lty = 1, col = 1,
             xlab = names[i], ylab = "Empirical CDF",
             main = "Weighted Empirical CDF")
        for (j in 1:(J - 1)) {
          w.ecdf.treat.j <- sapply(my.seq, weighted.ecdf,
                                   x = x$X[[i]][treat == j],
                                   weights = weights.mat[j + 1, treat == j])
          lines(my.seq, w.ecdf.treat.j, cex = 0.4,
                type = "s", pch = 16, lty = j + 1, col = j + 1)
        }
      }
    }
    par(ask = FALSE)
  }
}

