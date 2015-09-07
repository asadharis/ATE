\name{Cressie and Read objective functions}
\alias{Cressie and Read function}
\alias{cpp_cr_rho}
\alias{cpp_dcr_rho}
\alias{cpp_ddcr_rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cressie and Read objective functions
}
\description{
Cressie and Read class of objective functions parameterized by a scalar \eqn{\theta}.
}
\usage{
cpp_cr_rho(x,theta)
cpp_dcr_rho(x,theta)
cpp_ddcr_rho(x,theta)
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
\deqn{\rho(x) = -(1+\theta x)^{(\theta+1)/\theta}/(\theta+1) ,}
where \eqn{\theta} is scalar.  The expression is interpreted as limits for \eqn{\theta = 0} or  \eqn{\theta = -1} (see Newey and Smith (2004)).

For this package we use the modified version
\deqn{\rho(x) = -(1-\theta x)^{(\theta+1)/\theta}/(\theta+1) ,}
which led to better performance based on empirical simulations.
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
#Cressie-Read Family example
data("nsw")
Y<-nsw$re78[301:500]
X<-nsw[301:500,-c(1,9)]
treat<-nsw$treat[301:500]

#Exponential tilting
ATE(Y,treat,X,theta = 0)

#Empirical likelihood
ATE(Y,treat,X,theta = -1)

#Quadratic loss
ATE(Y,treat,X, theta = 1)
}
