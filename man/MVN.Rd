\name{MVN}
\alias{MVN}
\alias{dMVN}
\alias{rMVN}
\alias{rcMVN}
\title{
  Multivariate normal distribution
}
\description{
  Density and random generation for the multivariate normal distribution
  with mean equal to \code{mean}, precision matrix equal to \code{Q} (or covariance
  matrix equal to \code{Sigma}).

  Function \code{rcMVN} samples from the multivariate normal
  distribution with a canonical mean \eqn{b}, i.e., the mean is
  \eqn{\mu = Q^{-1}\,b.}{mu = Q^{-1} * b.}
}
\usage{
dMVN(x, mean=0, Q=1, Sigma, log=FALSE)

rMVN(n, mean=0, Q=1, Sigma)

rcMVN(n, b=0, Q=1, Sigma)
}
\arguments{
  \item{mean}{vector of mean.}
  \item{b}{vector of a canonical mean.}
  \item{Q}{precision matrix of the multivariate normal
    distribution. Ignored if \code{Sigma} is given.}
  \item{Sigma}{covariance matrix of the multivariate normal
    distribution. If \code{Sigma} is supplied, precision
    is computed from \eqn{\Sigma}{Sigma} as
    \eqn{Q = \Sigma^{-1}}{Q = Sigma^{-1}}.}
  \item{n}{number of observations to be sampled.}
  \item{x}{vector or matrix of the points where the density should be
    evaluated.}
  \item{log}{logical; if \code{TRUE}, log-density is computed}
}
\value{
Some objects.
}
\section{Value for dMVN}{
  A vector with evaluated values of the (log-)density
}  
\section{Value for rMVN}{
  A list with the components:
  \describe{
    \item{x}{vector or matrix with sampled values}
    \item{log.dens}{vector with the values of the log-density evaluated
      in the sampled values}
  }
}
\section{Value for rcMVN}{
  A list with the components:
  \describe{
    \item{x}{vector or matrix with sampled values}
    \item{mean}{vector or the mean of the normal distribution}    
    \item{log.dens}{vector with the values of the log-density evaluated
      in the sampled values}
  }
}
\references{
  Rue, H. and Held, L. (2005).
  \emph{Gaussian Markov Random Fields: Theory and Applications}.
  Boca Raton: Chapman and Hall/CRC.  
}
\seealso{
  \code{\link{dnorm}}, \code{\link[mvtnorm]{Mvnorm}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

### Univariate normal distribution
### ==============================
c(dMVN(0), dnorm(0))
c(dMVN(0, log=TRUE), dnorm(0, log=TRUE))

rbind(dMVN(c(-1, 0, 1)), dnorm(c(-1, 0, 1)))
rbind(dMVN(c(-1, 0, 1), log=TRUE), dnorm(c(-1, 0, 1), log=TRUE))

c(dMVN(1, mean=1.2, Q=0.5), dnorm(1, mean=1.2, sd=sqrt(2)))
c(dMVN(1, mean=1.2, Q=0.5, log=TRUE), dnorm(1, mean=1.2, sd=sqrt(2), log=TRUE))

rbind(dMVN(0:2, mean=1.2, Q=0.5), dnorm(0:2, mean=1.2, sd=sqrt(2)))
rbind(dMVN(0:2, mean=1.2, Q=0.5, log=TRUE), dnorm(0:2, mean=1.2, sd=sqrt(2), log=TRUE))

### Multivariate normal distribution
### ================================
mu <- c(0, 6, 8)
L <- matrix(1:9, nrow=3)
L[upper.tri(L, diag=FALSE)] <- 0
Sigma <- L \%*\% t(L)
Q <- chol2inv(chol(Sigma))
b <- solve(Sigma, mu)

dMVN(mu, mean=mu, Q=Q)
dMVN(mu, mean=mu, Sigma=Sigma)
dMVN(mu, mean=mu, Q=Q, log=TRUE)
dMVN(mu, mean=mu, Sigma=Sigma, log=TRUE)

xx <- matrix(c(0,6,8, 1,5,7, -0.5,5.5,8.5, 0.5,6.5,7.5), ncol=3, byrow=TRUE)
dMVN(xx, mean=mu, Q=Q)
dMVN(xx, mean=mu, Sigma=Sigma)
dMVN(xx, mean=mu, Q=Q, log=TRUE)
dMVN(xx, mean=mu, Sigma=Sigma, log=TRUE)

zz <- rMVN(1000, mean=mu, Sigma=Sigma)
rbind(apply(zz$x, 2, mean), mu)
var(zz$x)
Sigma
cbind(dMVN(zz$x, mean=mu, Sigma=Sigma, log=TRUE), zz$log.dens)[1:10,]

zz <- rcMVN(1000, b=b, Sigma=Sigma)
rbind(apply(zz$x, 2, mean), mu)
var(zz$x)
Sigma
cbind(dMVN(zz$x, mean=mu, Sigma=Sigma, log=TRUE), zz$log.dens)[1:10,]

zz <- rMVN(1000, mean=rep(0, 3), Sigma=Sigma)
rbind(apply(zz$x, 2, mean), rep(0, 3))
var(zz$x)
Sigma
cbind(dMVN(zz$x, mean=rep(0, 3), Sigma=Sigma, log=TRUE), zz$log.dens)[1:10,]


### The same using the package mvtnorm
### ==================================
# require(mvtnorm)
# c(dMVN(mu, mean=mu, Sigma=Sigma), dmvnorm(mu, mean=mu, sigma=Sigma))
# c(dMVN(mu, mean=mu, Sigma=Sigma, log=TRUE), dmvnorm(mu, mean=mu, sigma=Sigma, log=TRUE))
#
# rbind(dMVN(xx, mean=mu, Sigma=Sigma), dmvnorm(xx, mean=mu, sigma=Sigma))
# rbind(dMVN(xx, mean=mu, Sigma=Sigma, log=TRUE), dmvnorm(xx, mean=mu, sigma=Sigma, log=TRUE))
}
\keyword{distribution}
\keyword{multivariate}
