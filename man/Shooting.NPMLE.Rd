\name{expdata}
\alias{expdata}
\docType{data}
\title{
Example Data for Length-Biased Survival Models
}
\description{
An R package implementing the proposed estimation and variable selection procedures has been developed to accompany the paper titled ``Estimation and Variable Selection for Semiparametric Transformation Models with Length-Biased Survival Data.''

This dataset is provided to demonstrate the use of the package functions.
}
\usage{
data(expdata)
}
\format{
A data frame with \( n \) rows and multiple columns. The variables are:
\itemize{
  \item \code{A}: Left truncation time.
  \item \code{Y}: Observed survival time.
  \item \code{D}: Event indicator, 1 for failure and 0 for censored.
  \item \code{X}: Covariate matrix. In this example, we only use the first 10 columns for estimation.
}
}
\examples{
data(expdata)

A <- expdata[, 1]
Y <- expdata[, 2]
D <- expdata[, 3]
X <- expdata[, -c(1:3)]
X <- X[, 1:10]

ans0 <- Shooting.NPMLE(A, Y, D, X, 0.4, 10)
str(ans0)
}
\keywords{datasets}
