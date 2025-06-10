\name{Shooting.NPMLE}
\alias{Shooting.NPMLE}
\title{Estimation and Variable Selection for Semiparametric Transformation Models with Length-Biased Survival Data}
\usage{
Shooting.NPMLE(A, Y, D, X, r, nlan=10)
}
\arguments{
  \item{A}{Left truncation time.}
  \item{Y}{Observed survival time.}
  \item{D}{Event indicator: 1 for failure, 0 for censored.}
  \item{X}{Covariate matrix. In this example, only the first \code{p} columns are used for estimation.}
  \item{r}{Tuning parameter for the generalized log transformation model.}
}
\description{
This function implements the proposed estimation and variable selection procedure in the paper entitled
\emph{``Estimation and Variable Selection for Semiparametric Transformation Models with Length-Biased Survival Data''}.
It is part of an R package developed to accompany the methodology described in the paper.
}
\examples{
data(expdata)

A <- data[, 1]
Y <- data[, 2]
D <- data[, 3]
X <- data[, -c(1:3)]
ans0 <- Shooting.NPMLE(A, Y, D, X, r = 0.4, nlan = 10)
str(ans0)
}
\keywords{survival, semiparametric, transformation models, variable selection, length-biased}
