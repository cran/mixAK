\name{MVT}
\alias{MVT}
\alias{dMVT}
\alias{rMVT}
\alias{C_rMVT1_R}
\alias{C_dMVT1_R}
\alias{C_deriv_ldMVT_x}
\title{
  Multivariate Student t distribution
}
\description{
  Density and random generation for the multivariate Student t distribution
  with location equal to \code{mu}, precision matrix equal to \code{Q} (or scale
  matrix equal to \code{Sigma}).

  Mentioned functions implement the multivariate Student t distribution
  with a density given by
  \deqn{%
    p(\boldsymbol{z}) = 
      \frac{\Gamma\bigl(\frac{\nu+p}{2}\bigr)}{\Gamma\bigl(\frac{\nu}{2}\bigr)\,\nu^{\frac{p}{2}}\,\pi^{\frac{p}{2}}}\,
      \bigl|\Sigma\bigr|^{-\frac{1}{2}}\,
      \Bigl\{1 + \frac{(\boldsymbol{z} -
	\boldsymbol{\mu})'\Sigma^{-1}(\boldsymbol{z} -
	\boldsymbol{\mu})}{\nu}\Bigr\}^{-\frac{\nu+p}{2}},        
  }{%
    p(z) = (Gamma((nu+p)/2)/(Gamma(nu/2) * nu^(p/2) * pi^(p/2))) *
    |Sigma|^(-1/2) * (1 + ((z - mu)'*Sigma^(-1)*(z - mu)) / nu)^(-(nu+p)/2),
  }  
  where \eqn{p} is the dimension, \eqn{\nu > 0}{nu > 0} degrees of
  freedom, \eqn{\boldsymbol{\mu}}{mu} the location parameter and
  \eqn{\Sigma}{Sigma} the scale matrix.

  For \eqn{\nu > 1}{nu > 1}, the mean in equal to
  \eqn{\boldsymbol{\mu}}{mu},
  for \eqn{\nu > 2}{nu > 2}, the covariance matrix is equal to
  \eqn{\frac{\nu}{\nu - 2}\Sigma}{(nu / (nu - 2)) * Sigma}.        
}
\usage{
dMVT(x, df, mu=0, Q=1, Sigma, log=FALSE)

rMVT(n, df, mu=0, Q=1, Sigma)
}
\arguments{
  \item{df}{degrees of freedom of the multivariate Student t distribution.}
  \item{mu}{vector of the location parameter.}
  \item{Q}{precision (inverted scale) matrix of the multivariate Student
    t distribution. Ignored if \code{Sigma} is given.}
  \item{Sigma}{scale matrix of the multivariate Student t
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
\section{Value for dMVT}{
  A vector with evaluated values of the (log-)density
}  
\section{Value for rMVT}{
  A list with the components:
  \describe{
    \item{x}{vector or matrix with sampled values}
    \item{log.dens}{vector with the values of the log-density evaluated
      in the sampled values}
  }
}
%\references{
%}
\seealso{
  \code{\link{dt}}, \code{\link[mvtnorm]{Mvt}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

### Univariate central t distribution
z <- rMVT(10, df=1, mu=0, Q=1)
ldz <- dMVT(z$x, df=1, log=TRUE)
boxplot(as.numeric(z$x))
cbind(z$log.dens, ldz, dt(as.numeric(z$x), df=1, log=TRUE))

### Multivariate t distribution
mu <- c(1, 2, 3)
Sigma <- matrix(c(1, 1, -1.5,  1, 4, 1.8,  -1.5, 1.8, 9), nrow=3)
Q <- chol2inv(chol(Sigma))

nu <- 3
z <- rMVT(1000, df=nu, mu=mu, Sigma=Sigma)
apply(z$x, 2, mean)              ## should be close to mu
((nu - 2) / nu) * var(z$x)       ## should be close to Sigma            

dz <- dMVT(z$x, df=nu, mu=mu, Sigma=Sigma)
ldz <- dMVT(z$x, df=nu, mu=mu, Sigma=Sigma, log=TRUE)

### Compare with mvtnorm package
#require(mvtnorm)
#ldz2 <- dmvt(z$x, sigma=Sigma, df=nu, delta=mu, log=TRUE)
#plot(z$log.dens, ldz2, pch=21, col="red3", bg="orange", xlab="mixAK", ylab="mvtnorm")
#plot(ldz, ldz2, pch=21, col="red3", bg="orange", xlab="mixAK", ylab="mvtnorm")
}
\keyword{distribution}
\keyword{multivariate}
