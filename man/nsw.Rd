\name{nsw}
\alias{nsw}
\docType{data}
\title{
Lalonde Subsample of the National Supported Work Demonstration Data (NSW)
}
\description{
This data was collected to evaluate the National Supported Work (NSW) Demonstration project in Lalonde (1986).
}
\usage{data(nsw)}

\format{
  A data frame with 722 observations on the following 9 variables.
  \describe{
    \item{\code{treat}}{Vector indicating treatment assignment}
    \item{\code{age}}{Age of subject}
    \item{\code{ed}}{Years of schooling}
    \item{\code{black}}{Black race (1 if Black, 0 otherwise)}
    \item{\code{hisp}}{Hispanic race (1 if Hispanic, 0 otherwise)}
    \item{\code{married}}{Marital status (1 if married, 0 otherwise)}
    \item{\code{nodeg}}{College degree (1 if no college degree, 0 otherwise)}
    \item{\code{re75}}{Earnings in 1975}
    \item{\code{re78}}{The response variable. Earnings in 1978}
  }
}
\details{
This data set was obtained from "NSW Data Files" and is a subset of the data collected for the original National Supported Work Evaluation Study (see source below). A portion of the summary of the study is quoted below. This quote was obtained from the website of the "Interuniversity Consortium for Political and Social Research" (see source below).

"This study is an evaluation of the National Supported Work Demonstration project, a transitional, subsidized work experience program for four target groups of people with longstanding employment problems: ex-offenders, former drug addicts, women who were long-term recipients of welfare benefits, and school dropouts, many with criminal records. The program provided up to 12-18 months of employment to about 10,000 individuals at 15 locations across the country for four years. In ten of these sites -- Atlanta, Chicago, Hartford, Jersey City, Newark, New York, Philadelphia, Oakland, San Francisco, and Wisconsin, 6,600 eligible applicants were randomly assigned either to experimental groups (offered a job in supported work) or to control groups, and an evaluation was conducted on the effects of the Supported Work Program. "
}
\source{
"NSW Data Files" from Rajeev Dehejia's website. URL: \url{http://users.nber.org/~rdehejia/data/nswdata2.html}

"National Supported Work Evaluation Study, 1975-1979: Public Use Files." from the Interuniversity Consortium for Political and Social Research. URL: \url{http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/7865}

}

\references{

Lalonde, R. (1986). "Evaluating the Econometric Evaluations of Training Programs," American Economic Review, 76(4), 604-620.

Dehejia R. and Wahba S. (1999). "Causal Effects in Non-Experimental Studies: Re-Evaluating the Evaluation of Training Programs," JASA 94 (448), 1053-1062.

}


\examples{
library(ATE)
# NSW Data example
data("nsw")

Y <- nsw$re78
X <- nsw[, c(-1, -9)]
treat <- nsw$treat

fit1 <- ATE(Y, treat, as.matrix(X))
fit1
summary(fit1)
# plot(fit1)
}
\keyword{datasets}
