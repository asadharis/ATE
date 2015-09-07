\name{summary.ATE}
\alias{summary.ATE}
\alias{print.summary.ATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing output of study.
}
\description{
\code{summary} method for class \code{"ATE"}
}
\usage{
\method{summary}{ATE}(object, ...)

\method{print}{summary.ATE}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{"ATE"}, usually a result of a call 
to \code{\link{ATE}}.
}
\item{x}{
An object of class \code{"summary.ATE"}, usually a result of 
a call to \code{\link{summary.ATE}}.
}
  \item{\dots}{
Further arguments passed to or from methods.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{print.summary.ATE} prints a simplified output similar to \code{\link{print.summary.lm}}. The resulting table provides the point estimates, estimated standard errors, 95\% Wald confidence intervals, the Z-statistic and the P-values for a Z-test.
}
\value{
The function \code{summary.ATE} returns a list with the following components
\item{Estimate}{
A matrix with point estimates along with standard errors, confidence intervals etc. This is the matrix users see with the \code{print.summary.RIPW} function.
}
\item{vcov}{
The variance-covariance matrix of the point estimates.
}
\item{Conv}{
The convergence result of the \code{object}.
}
\item{weights}{
The weights for each subject in each treatment arm. These are same as the \code{weight} component of the \code{"RIPW" object}.
}
\item{call}{
The \code{call} passed on as an argument of the function which is equivalent to \code{object$call}.
}
}
\author{
Asad Haris, Gary Chan
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ATE}}
}
\examples{
library(ATE)
#binary treatment

set.seed(25)
n <- 200
Z <- matrix(rnorm(4*n),ncol=4,nrow=n)
prop <- 1 / (1 + exp(Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))
treat <- rbinom(n, 1, prop)
Y <- 200 + 10*treat+ (1.5*treat-0.5)*(27.4*Z[,1] + 13.7*Z[,2] +
          13.7*Z[,3] + 13.7*Z[,4]) + rnorm(n)
X <- cbind(exp(Z[,1])/2,Z[,2]/(1+exp(Z[,1])),
          (Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)

#estimation of average treatment effects (ATE)
fit1<-ATE(Y,treat,X)
summary(fit1)
#plot(fit1)

#estimation of average treatment effects on treated (ATT)
fit2<-ATE(Y,treat,X,ATT=TRUE)
summary(fit2)
#plot(fit2)

#three treatment groups
set.seed(25)
n <- 200
Z <- matrix(rnorm(4*n),ncol=4,nrow=n)
prop1 <- 1 / (1 + exp(1+Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))
prop2 <- 1 / (1 + exp(Z[,1] - 0.5 * Z[,2] + 0.25*Z[,3] + 0.1 * Z[,4]))

U <-runif(n)
treat <- numeric(n)
treat[U>(1-prop2)]=2
treat[U<(1-prop2)& U>(prop2-prop1)]=1

Y <- 210 + 10*treat +(27.4*Z[,1] + 13.7*Z[,2] + 
            13.7*Z[,3] + 13.7*Z[,4]) + rnorm(n)
X <- cbind(exp(Z[,1])/2,Z[,2]/(1+exp(Z[,1])),
            (Z[,1]*Z[,3]/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)

fit3<-ATE(Y,treat,X)
summary(fit3)
#plot(fit3)

}
