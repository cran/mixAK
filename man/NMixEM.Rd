\name{NMixEM}
\alias{NMixEM}
\alias{print.NMixEM}
\title{
  EM algorithm for a homoscedastic normal mixture
}
\description{
  This function computes ML estimates of the parameters
  of the \eqn{p}-dimensional \eqn{K}-component normal mixture using the EM algorithm
}
\usage{
NMixEM(y, K, weight, mean, Sigma, toler=1e-5, maxiter=500)

\method{print}{NMixEM}(x, \dots)
}
\arguments{
  \item{y}{vector (if \eqn{p=1}) matrix or data frame (if \eqn{p > 1}) with data. Rows correspond to
    observations, columns correspond to margins.}
  \item{K}{required number of mixture components.}
  \item{weight}{a numeric vector with initial mixture weights.

    If not given, initial weights are all equal to \eqn{1/K}.
  }
  \item{mean}{vector or matrix of initial mixture means.
    
    For \eqn{p=1} this
    should be a vector of length \eqn{K}, for \eqn{p>1} this should be a
    \eqn{K\times p}{K x p} matrix with mixture means in rows.
  }
  \item{Sigma}{number or \eqn{p\times p}{p x p} matrix giving the
    initial variance/covariance matrix.
  }
  \item{toler}{tolerance to determine convergence.}
  \item{maxiter}{maximum number of iterations of the EM algorithm.}
  \item{x}{an object of class \code{NMixEM}.}
  \item{\dots}{additional arguments passed to the default \code{print} method.}  
}
\value{
  An object of class \code{NMixEM} which has the following components:
  \item{K}{number of mixture components}
  \item{weight}{estimated mixture weights}
  \item{mean}{estimated mixture means}
  \item{Sigma}{estimated covariance matrix}
  \item{loglik}{log-likelihood value at fitted values}
  \item{aic}{Akaike information criterion
    (\eqn{-2\hat{\ell} + 2\nu}{-2*loglik + 2*nu}), where
    \eqn{\hat{\ell}}{loglik} stands for the log-likelihood value at
    fitted values and \eqn{\nu}{nu} for the number of free model parameters
  }
  \item{bic}{Bayesian (Schwarz) information criterion
    (\eqn{-2\hat{\ell} + \log(n)\,\nu}{-2*loglik + log(n)*nu}), where
    \eqn{\hat{\ell}}{loglik} stands for the log-likelihood value at
    fitted values and \eqn{\nu}{nu} for the number of free model
    parameters, and \eqn{n} for the sample size    
  }  
  \item{iter}{number of iterations of the EM algorithm used to get the solution}
  \item{iter.loglik}{values of the log-likelihood at iterations of the
    EM algorithm}
  \item{iter.Qfun}{values of the EM objective function at iterations of
    the EM algorithm}
  \item{dim}{dimension \eqn{p}}
  \item{nobs}{number of observations \eqn{n}}       
}
\references{
  Dempster, A. P., Laird, N. M., Rubin, D. B. (1977).
  Maximum likelihood from incomplete data via the EM algorithm.
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{39}, 1-38.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
\dontrun{
## Estimates for 3-component mixture in  Anderson's iris data
## ==========================================================
data(iris, package="datasets")
summary(iris)

VARS <- names(iris)[1:4]
fit <- NMixEM(iris[, VARS], K = 3)
print(fit)

apply(subset(iris, Species == "versicolor")[, VARS], 2, mean)
apply(subset(iris, Species == "setosa")[, VARS], 2, mean)
apply(subset(iris, Species == "virginica")[, VARS], 2, mean)

## Estimates of 6-component mixture in Galaxy data
## ==================================================
data(Galaxy, package="mixAK")
summary(Galaxy)

fit2 <- NMixEM(Galaxy, K = 6)
y <- seq(5, 40, length=300)
fy <- dMVNmixture(y, weight=fit2$weight, mean=fit2$mean,
                     Sigma=rep(fit2$Sigma, fit2$K))
hist(Galaxy, prob=TRUE, breaks=seq(5, 40, by=0.5),
     main="", xlab="Velocity (km/sec)", col="sandybrown")
lines(y, fy, col="darkblue", lwd=2)
}
}
\keyword{multivariate}
