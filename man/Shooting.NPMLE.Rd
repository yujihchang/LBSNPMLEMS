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
\Outputs{
\item{npmle}{Nonparametric maximum likelihood estimator without variable selection.}
\item{haz}{Estimated cumulative baseline hazard function.}
\item{npmle.sd}{Estimated standard deviation of the nonparametric maximum likelihood estimator.}
\item{lambda.aic}{Selected tuning parameter \code{lambda} based on AIC.}
\item{lambda.bic}{Selected tuning parameter \code{lambda} based on BIC.}
\item{AIC}{Akaike information criterion value.}
\item{BIC}{Bayesian information criterion value.}
\item{coef.aic}{Regression coefficient estimates selected based on AIC.}
\item{coef.bic}{Regression coefficient estimates selected based on BIC.}
\item{coef.aic.sd}{Estimated standard deviation of \code{coef.aic}.}
\item{coef.bic.sd}{Estimated standard deviation of \code{coef.bic}.}
\item{choose.X.aic}{Indices of covariates selected based on AIC.}
\item{choose.X.bic}{Indices of covariates selected based on BIC.}
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
