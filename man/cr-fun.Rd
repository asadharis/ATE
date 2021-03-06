\name{CRFamily}
\alias{Cressie and Read function}
\alias{CRFamily}
\alias{CRFamilyDerivative}
\alias{CRFamilySecondDerivative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cressie and Read Family of Functions
}
\description{
Cressie and Read class of functions parameterized by a scalar \eqn{\theta} used
within this package. In our case, we use a slight modification of the original
CR Family of functions, see details.
}
\usage{
CRFamily(x, theta)
CRFamilyDerivative(x, theta)
CRFamilySecondDerivative(x, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A vector of points at which the function will be evaluated at.
}
\item{theta}{
A user defined value for \eqn{\theta}{\theta}. This must be a scalar. See details.
}
}
\details{
This set of functions can be used for the main \code{\link{ATE}} function for the arguments \code{theta}. The class of functions, are defined by (see Cressie and Read (1984))
\deqn{\rho(x) = -\frac{ (1 + \theta x) ^ {(\theta + 1) / \theta} }{
        \theta+1} ,}{\rho(x) = -(1 + \theta x)^{(\theta + 1) / \theta}/(\theta + 1) ,}
where \eqn{\theta} is scalar.  The expression is interpreted as limits for \eqn{\theta = 0} or  \eqn{\theta = -1} (see Newey and Smith (2004)).

For this package we use the modified version
\deqn{\rho(x) = -\frac{(1-\theta x)^{(\theta+1)/\theta} }{\theta+1} , }{\rho(x) = -(1-\theta x)^{(\theta+1)/\theta}/(\theta+1) ,}
which leads to greater computational stability (based on empirical simulations).
}
\value{
The output of the real values function for a given parameter \eqn{\theta}{\theta}.
}
\references{
Chan, K. C. G. and Yam, S. C. P. and Zhang, Z. (2015). "Globally Efficient Nonparametric Infernence of Average Treatment Effects by Emperical Balancing Calibration Weighting", under revision.

Chan, K. C. G. and Yam, S. C. P. (2014). "Oracle, Multiple Robust and Multipurpose Calibration in a Missing Response Problem".  Statistical Science 29 (3): 380-386.

Cressie, N. and Read, T. R. C. (1984). "Multinomial goodness-of-fit tests".  J. Roy. Statist. Soc. Ser. B 46 440-464.

Newey, W. K. and Smith, R. J. (2004). Higher order properties
of GMM and generalized empirical likelihood estimators.
Econometrica 72 219-255.
}
\author{
Asad Haris, Gary Chan
}

\examples{
# Cressie-Read Family example
data("nsw")
Y <- nsw$re78
X <- nsw[, -c(1, 9)]
X <- as.matrix(X)
treat <- nsw$treat

# Exponential tilting.
ATE(Y, treat, X, theta = 0)

# Empirical likelihood.
ATE(Y, treat, X, theta = -1)

# Quadratic loss.
ATE(Y, treat, X, theta = 1)
}
