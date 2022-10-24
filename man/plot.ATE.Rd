\name{plot.ATE}
\alias{plot.ATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of Empirical CDF and Weighted Empirical CDF for Covariates
}
\description{
Plot function for class \code{'ATE'}.
}
\usage{
\method{plot}{ATE}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of type \code{'ATE'}.
}
  \item{\dots}{
Further arguments passed to or from the function.
}
}
\details{
This function plots the empirical CDF and weighted empirical CDF for each covariate to demonstrate the effect of covariate balancing and for graphical diagnostics. In  observational studies with confounding, the covariate distributions are different for each treatment arms.  Comparisons of unweighted empirical CDF would demonstrate this difference. The balancing weights constructed by balancing moments of covariate distributions, and the weighted CDF would show an improved balance.
}

\author{
Asad Haris, Gary Chan
}

\seealso{
\code{\link{ATE}}
}
\examples{
# Binary treatment.
set.seed(25)
n <- 200
Z <- matrix(rnorm(4 * n), ncol = 4, nrow = n)
prop <- 1 / (1 + exp(Z[, 1] - 0.5 * Z[, 2] + 0.25 * Z[, 3] + 0.1 * Z[, 4]))
treat <- rbinom(n, 1, prop)
Y <- 200 + 10 * treat + (1.5 * treat - 0.5) *
     (27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] + 13.7 * Z[, 4]) + rnorm(n)
X <- cbind(exp(Z[, 1]) / 2, Z[, 2] / (1 + exp(Z[, 1])),
           (Z[, 1] * Z[, 3] / 25 + 0.6) ^ 3, (Z[, 2] + Z[, 4] + 20) ^ 2)

# Estimation of average treatment effects (ATE).
fit1 <- ATE(Y, treat, X)
summary(fit1)
\dontrun{
plot(fit1)
}
# Estimation of average treatment effects on treated (ATT).
fit2 <- ATE(Y, treat, X, ATT = TRUE)
summary(fit2)
\dontrun{
plot(fit2)
}

# Four treatment groups.
set.seed(25)
n <- 200
Z <- matrix(rnorm(4 * n), ncol = 4, nrow = n)
prop1 <- 1 / (1 + exp(1 + Z[, 1] - 0.5 * Z[, 2] + 0.25 * Z[, 3] + 0.1 * Z[, 4]))
prop2 <- 1 / (1 + exp(Z[, 1] - 0.5 * Z[, 2] + 0.25 * Z[, 3] + 0.1 * Z[, 4]))

U <- runif(n)
treat <- numeric(n)
treat[U > (1 - prop2)] = 2
treat[U < (1 - prop2) & U > (prop2 - prop1)] = 1

Y <- 210 + 10 * treat + (27.4 * Z[, 1] + 13.7 * Z[, 2] + 13.7 * Z[, 3] +
                        13.7 * Z[, 4]) + rnorm(n)
X <- cbind(exp(Z[, 1]) / 2, Z[, 2] / (1 + exp(Z[, 1])),
           (Z[, 1] * Z[, 3] / 25 + 0.6) ^ 3, (Z[, 2] + Z[, 4] + 20) ^ 2)

fit3 <- ATE(Y, treat, X)
summary(fit3)
\dontrun{
plot(fit3)
}
}
