\name{ATE}
\alias{ATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the Average Treatment Effect
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The main function for estimating the average treatment effect or the average
treatment effect on the treated. This function creates an \code{ATE} object
which can be used as inputs for generic S3 \code{plot} or \code{summary}
functions. This function uses a covariate balancing method which creates
weights for each subject, without a need to specific a propensity score
or an outcome regression model. The main function implements the BFGS
algorithm with backtracking. For details of the estimator
see Chan et al. (2015).
}
\usage{
ATE(Y, treat, X, theta = 1, ATT = FALSE,
    verbose = FALSE, max.iter = 100, tol = 1e-10,
    initial.values = NULL,
    backtrack.alpha = 0.3,
    backtrack.beta = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{Y}{
The response vector of length \eqn{n}. This has to be a numeric vector.
}

 \item{treat}{
The vector of treatment assignments of length \eqn{n}. Must be
coded as \eqn{0, 1, \ldots}{\{0, 1, ...\}} and must be binary
vector when \code{ATT = TRUE}.
}
 \item{X}{
A \eqn{N \times p}{N x p}-matrix or data frame of covariates \code{X} to be balanced.
This matrix does not need to include an intercept.
}
\item{theta}{
A real scalar parameter for the Cressie-Read family of objective functions.
The default is \eqn{\theta = 0} (exponential tilting).  Other popular
examples are \eqn{\theta = -1} (empirical likelihood) and
\eqn{\theta = 1} (quadratic loss). See \code{\link{CRFamily}} for details.
}
\item{ATT}{
A logical value indicating whether to estimate the
treatment effect on the treated.
}
\item{verbose}{
A logical value indicating whether to print the progress of the function.
Default value is \code{FALSE}.
}
\item{max.iter}{
The maximum number of iterations for the BFGS algorithm.
}
\item{tol}{
The absolute tolerance used to determine a stopping criteria for
the BFGS algorithm.
}
\item{initial.values}{
A numeric vector or matrix of possible initial values for the BFGS algorithm.
Must be a \eqn{J\times K}{J x K} matrix where \eqn{J}{J} is the number
of treatment arms. For \code{ATT = TRUE} this must be a \eqn{K}{K}-vector.
}
\item{backtrack.alpha}{
A scalar parameter for the backtracking line search algorithm
with \eqn{\alpha \in (0,0.5)}{\alpha in (0,0.5)}.
}
\item{backtrack.beta}{
A scalar parameter for for the backtracking line search algorithm
with \eqn{\beta \in (0,1)}{\beta in (0,1)}.
}

}


\value{
%%  ~Describe the value returned
The function returns an object of type \code{'ATE'}, a list with the following elements:

\item{estimate}{
The vector of point estimates for the average treatment effect.
For a binary treatment it also contains the average difference of treatment effects.
}
\item{vcov}{
The estimated variance covariance matrix for the estimates of the
treatment effects for each treatment group.
}
\item{lambda}{
The resulting solution of the main optimization problems,
\eqn{\hat{\lambda}}{hat{\lambda}}, as described in Chan et al.(2015).
In the case of a simple, binary treatment study, the object has
\code{lambda.treat} and \code{lambda.placebo} and when \code{ATT = TRUE},
we only have \code{lambda.placebo}. For a multiple treatment study design
we have \code{lam.mat}, a matrix with each row for each
\eqn{\hat{\lambda}}{hat{\lambda}}, corresponding to each treatment arm.
}
\item{weights}{
The weights obtained by the balancing covariate method for each
treatment group. In the case of \code{ATT = TRUE}, we only have weights for
the untreated. For binary treatment the list will have \code{weights.treat} and
\code{weights.placebo}. For multiple treatment effect the list contains
a \eqn{J\times N} matrix \code{weights.mat}.
}
\item{gp}{
A string specifying the type of study design. For binary treatment effect
with \code{ATT = FALSE} is denoted by group "simple".
With \code{ATT = TRUE} we have "ATT" and finally "MT" is for
multiple treatment arms.
}
\item{converge}{
A logical value indicating convergence of BFGS algorithm.
}
\item{X, Y, treat}{
The data used for estimation.
}
\item{J}{
A scalar indicating the number of treatment arms.
}
\item{K}{
A scalar indicating the one plus the dimension of the range space of X.
\eqn{K = 1 + p}.
}
\item{call}{
The matched call.
}

}

\references{
Chan, K. C. G., Yam, S. C. P., & Zhang, Z. (2016). Globally efficient non-parametric inference of average treatment effects by empirical balancing calibration weighting. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(3), 673-700. \doi{10.1111/rssb.12129}

Haris, Asad and Chan, K. C. G. (2015). "\pkg{ATE}: An \code{R} package for Nonparametric Inference of Average Treatment Effects", under revision
}

\author{
Asad Haris, Gary Chan.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.ATE}},
\code{\link{plot.ATE}}
}
\examples{
library(ATE)
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
# plot(fit1)

# Estimation of average treatment effects on treated (ATT).
fit2 <- ATE(Y, treat, X, ATT = TRUE)
summary(fit2)
# plot(fit2)


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
# plot(fit3)
}
