\name{nhanes_bmi}
\alias{nhanes_bmi}
\docType{data}
\title{
Subsample of National Health and Nutritional Examination Survey for School Children
}
\description{
This data is a limited subset of the 2007-2008 National Health and Nutrition Examination Survey (NHANES), for studying whether particpation of school meal programs lead to an increase in body mass index (BMI).
}
\usage{data(nhanes_bmi)}

\format{
  A data frame with 2330 observations on the following 13 variables.
  \describe{
    \item{\code{BMI}}{The response variable.  Body mass index.}
    \item{\code{School_meal}}{Participation of school meal programs (1: Yes, 0: No)}
    \item{\code{age}}{Age of Child}
    \item{\code{ChildSex}}{Gender of the Child (1: Male, 0: Female)}
    \item{\code{black}}{Black Race (1: Black, 0: otherwise)}
    \item{\code{mexam}}{Hispanic Race (1: Hispanic: 0 otherwise)}
    \item{\code{pir200_plus}}{Family above 200\% of the federal poverty level  (1: Yes, 0: No)}
    \item{\code{WIC}}{Participation in Special Supplemental Nutrition program (1: Yes, 0: No)}
    \item{\code{Food_Stamp}}{Participation in foold stamp program (1: Yes, 0: No)}
    \item{\code{fsdchbi}}{Childhood food security (1: Secure, 0: insecure)}
    \item{\code{AnyIns}}{Any insurance (1: Yes, 0: No)}
    \item{\code{RefSex}}{Gender of the adult respondent (1: Male, 0: Female)}
    \item{\code{RefAge}}{Age of the adult respondent}
  }
}

\source{
Data adapted from the 2007-2008 National Health and Nutrition Examination Survey (NHANES).
\url{http://www.cdc.gov/nchs/nhanes.htm}
}


\examples{
library(ATE)
data("nhanes_bmi")

Y <- nhanes_bmi[, 1]
treat <- nhanes_bmi[, 2]
X <- nhanes_bmi[, -c(1, 2)]
X <- as.matrix(X)

fit1 <- ATE(Y, treat, X)
fit1
summary(fit1)
# plot(fit1)
}
\keyword{datasets}
