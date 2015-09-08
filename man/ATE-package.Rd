\name{ATE-package}
\alias{ATE-package}
\docType{package}
\title{
Nonparametric efficient inference of the average treatment effects for observational data
}
\description{
This package provides a user-friendly interface for nonparametric efficient inference of average treatment effects for observational data. The package provides point estimates for average treatment effects, average treatment effect on the treated and can also handle the case of multiple treatments. The package also allows inference by consistent variance estimates. The point estimates for the treatment effect and variance estimates are described in Chan et al. (2015).
}
\details{

\tabular{ll}{
Package: \tab ATE\cr
Type: \tab Package\cr
Version: \tab 0.3.2\cr
Date: \tab 2015-09-08\cr
License: \tab GPL (>= 2) \cr
}

The package includes the following functions:

\tabular{ll}{
\code{\link{ATE}}: \tab Estimate the average treatment effect \cr
\code{\link{summary.ATE}}: \tab \code{summary} method for class \code{"ATE"} \cr
\code{\link{plot.ATE}}: \tab Plot function for class "RIPW" \cr
}
}

\author{
Asad Haris, Gary Chan

Maintainer: Asad Haris <aharis@uw.edu>
}

\references{
Chan, K. C. G. and Yam, S. C. P. and Zhang, Z. (2015). "Globally Efficient Nonparametric Inference of Average Treatment Effects by Empirical Balancing Calibration Weighting", under review.
}

\keyword{ package }
\seealso{
\code{\link{ATE}},
\code{\link{summary.ATE}},
\code{\link{plot.ATE}},
\code{\link{cpp_cr_rho}}
}
