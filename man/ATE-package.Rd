\name{ATE-package}
\alias{ATE-package}
\docType{package}
\title{
Nonparametric efficient inference of the average treatment effects for observational data
}
\description{
This package provides a user-friendly interface for nonparametric efficient inference of average treatment effects for observational data. The package provides point estimates for average treatment effects, average treatment effect on the treated and can also handle the case of multiple treatments. The package also allows inference by consistent variance estimates. The point estimates for the treatment effect and variance estimates are described in Chan et al. (2016).
}
\details{
The package includes the following functions:

\tabular{ll}{
\code{\link{ATE}}: \tab Estimate the average treatment effect \cr
\code{\link{summary.ATE}}: \tab \code{summary} method for class \code{'ATE'} \cr
\code{\link{plot.ATE}}: \tab Plot function for class \code{'ATE'} \cr
}
}

\author{
Asad Haris, Gary Chan

Maintainer: Asad Haris <aharis@uw.edu>
}

\references{
Chan, K. C. G., Yam, S. C. P., & Zhang, Z. (2016). Globally efficient non-parametric inference of average treatment effects by empirical balancing calibration weighting. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(3), 673-700. \doi{10.1111/rssb.12129}
}

\keyword{ package }
\seealso{
\code{\link{ATE}},
\code{\link{summary.ATE}},
\code{\link{plot.ATE}},
\code{\link{CRFamily}}
}
