\name{MVNmixture}
\alias{MVNmixture}
\alias{dMVNmixture}
\alias{dMVNmixture2}
\alias{rMVNmixture}
\alias{rMVNmixture2}
\title{
  Mixture of (multivariate) normal distributions
}
\description{
  Density and random generation for the mixture of the \eqn{p}-variate normal distributions
  with means given by \code{mean}, precision matrix given by \code{Q} (or covariance
  matrices given by\code{Sigma}).  
}
\usage{
dMVNmixture(x, weight, mean, Q, Sigma, log=FALSE)

dMVNmixture2(x, weight, mean, Q, Sigma, log=FALSE)

rMVNmixture(n, weight, mean, Q, Sigma)

rMVNmixture2(n, weight, mean, Q, Sigma)
}
\details{
  Functions \code{dMVNmixture} and \code{dMVNmixture2} differ only
  internally in the way they compute the mixture density. In
  \code{dMVNmixture}, only multivariate normal densities are evaluated
  in compiled C++ code and mixing is done directly in R. In
  \code{dMVNmixture2}, everything is evaluated in compiled C++
  code. Normally, both \code{dMVNmixture} and \code{dMVNmixture2} should
  return the same results.

  Similarly for \code{rMVNmixture} and \code{rMVNmixture2}. Another
  difference is that \code{rMVNmixture} returns only random generated
  points and \code{rMVNmixture2} also values of the density evaluated in
  the generated points.
}  
\arguments{
  \item{weight}{vector of length \eqn{K} with the mixture weights or
    values which are proportional to the weights.}
  \item{mean}{vector or matrix of mixture means.

    For \eqn{p=1} this
    should be a vector of length \eqn{K}, for \eqn{p>1} this should be a
    \eqn{K\times p}{K x p} matrix with mixture means in rows.
  }
  \item{Q}{precision matrices of the multivariate normal
    distribution. Ignored if \code{Sigma} is given.

    For \eqn{p=1} this should be a vector of length \eqn{K}, for \eqn{p
      > 1} this should be a list of length \eqn{K} with the mixture
    precision matrices as components of the list.
  }
  \item{Sigma}{covariance matrix of the multivariate normal
    distribution. If \code{Sigma} is supplied, precisions
    are computed from \eqn{\Sigma}{Sigma} as
    \eqn{Q = \Sigma^{-1}}{Q = Sigma^{-1}}.

    For \eqn{p=1} this should be a vector of length \eqn{K}, for \eqn{p
      > 1} this should be a list of length \eqn{K} with the mixture
    covariance matrices as components of the list.
  }
  \item{n}{number of observations to be sampled.}
  \item{x}{vector or matrix of the points where the density should be
    evaluated.}
  \item{log}{logical; if \code{TRUE}, log-density is computed}
}
\value{
Some objects.
}
\section{Value for dMVNmixture}{
  A vector with evaluated values of the (log-)density.
}  
\section{Value for dMVNmixture2}{
  A vector with evaluated values of the (log-)density.
}  
\section{Value for rMVNmixture}{
  A vector (for \code{n=1} or for univariate mixture)
  or matrix with sampled values (in rows of the matrix).
}
\section{Value for rMVNmixture2}{
  A list with components named \code{x} which is 
  a vector (for \code{n=1} or for univariate mixture)
  or matrix with sampled values (in rows of the matrix)
  and \code{dens} which are the values of the density evaluated in \code{x}.
}
\seealso{
  \code{\link{dnorm}}, \code{\link{MVN}}, \code{\link[mvtnorm]{Mvnorm}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

##### Univariate normal mixture
##### =========================
mu <- c(-1, 1)
Sigma <- c(0.25^2, 0.4^2)
Q <- 1/Sigma
w <- c(0.3, 0.7)

xx <- seq(-2, 2.5, length=100)
yyA <- dMVNmixture(xx, weight=w, mean=mu, Sigma=Sigma)
yyB <- dMVNmixture(xx, weight=w, mean=mu, Q=Q)
yyC <- dMVNmixture2(xx, weight=w, mean=mu, Sigma=Sigma)
yyD <- dMVNmixture2(xx, weight=w, mean=mu, Q=Q)

xxSample <- rMVNmixture(1000, weight=w, mean=mu, Sigma=Sigma)
xxSample2 <- rMVNmixture2(1000, weight=w, mean=mu, Sigma=Sigma)

sum(abs(xxSample2$dens - dMVNmixture(xxSample2$x, weight=w, mean=mu, Sigma=Sigma)) > 1e-15)
xxSample2 <- xxSample2$x

par(mfrow=c(2, 2), bty="n")
plot(xx, yyA, type="l", col="red", xlab="x", ylab="f(x)")
points(xx, yyB, col="darkblue")
hist(xxSample, col="lightblue", prob=TRUE, xlab="x", xlim=range(xx), ylim=c(0, max(yyA)),
     main="Sampled values")
lines(xx, yyA, col="red")
plot(xx, yyC, type="l", col="red", xlab="x", ylab="f(x)")
points(xx, yyD, col="darkblue")
hist(xxSample2, col="sandybrown", prob=TRUE, xlab="x", xlim=range(xx), ylim=c(0, max(yyA)),
     main="Sampled values")
lines(xx, yyC, col="red")


##### Bivariate normal mixture
##### ========================
### Choice 1
sd11 <- sd12 <- 1.1
sd21 <- 0.5
sd22 <- 1.5
rho2 <- 0.7
Xlim <- c(-3, 4)
Ylim <- c(-6, 4)

### Choice 2
sd11 <- sd12 <- 0.3
sd21 <- 0.5
sd22 <- 0.3
rho2 <- 0.8
Xlim <- c(-3, 2.5)
Ylim <- c(-2.5, 2.5)

mu <- matrix(c(1,1, -1,-1), nrow=2, byrow=TRUE)
Sigma <- list(diag(c(sd11^2, sd12^2)),
              matrix(c(sd21^2, rho2*sd21*sd22, rho2*sd21*sd22, sd22^2), nrow=2))
Q <- list(chol2inv(chol(Sigma[[1]])), chol2inv(chol(Sigma[[2]])))
w <- c(0.3, 0.7)

xx1 <- seq(mu[2,1]-3*sd21, mu[1,1]+3*sd11, length=100)
xx2 <- seq(mu[2,2]-3*sd22, mu[1,2]+3*sd12, length=90)
XX <- cbind(rep(xx1, length(xx2)), rep(xx2, each=length(xx1)))
yyA <- matrix(dMVNmixture(XX, weight=w, mean=mu, Sigma=Sigma), nrow=length(xx1), ncol=length(xx2))
yyB <- matrix(dMVNmixture(XX, weight=w, mean=mu, Q=Q), nrow=length(xx1), ncol=length(xx2))
yyC <- matrix(dMVNmixture2(XX, weight=w, mean=mu, Sigma=Sigma), nrow=length(xx1), ncol=length(xx2))
yyD <- matrix(dMVNmixture2(XX, weight=w, mean=mu, Q=Q), nrow=length(xx1), ncol=length(xx2))

#xxSample <- rMVNmixture(1000, weight=w, mean=mu, Sigma=Sigma)
### Starting from version 3.6, the above command led to SegFault
### on CRAN r-patched-solaris-sparc check.
### Commented here on 20140806 (version 3.6-1).
xxSample2 <- rMVNmixture2(1000, weight=w, mean=mu, Sigma=Sigma)

sum(abs(xxSample2$dens - dMVNmixture(xxSample2$x, weight=w, mean=mu, Sigma=Sigma)) > 1e-15)
xxSample2 <- xxSample2$x

par(mfrow=c(1, 2), bty="n")
plot(xxSample, col="darkblue", xlab="x1", ylab="x2", xlim=Xlim, ylim=Ylim)
contour(xx1, xx2, yyA, col="red", add=TRUE)
plot(xxSample2, col="darkblue", xlab="x1", ylab="x2", xlim=Xlim, ylim=Ylim)
contour(xx1, xx2, yyC, col="red", add=TRUE)

par(mfrow=c(2, 2), bty="n")
contour(xx1, xx2, yyA, col="red", xlab="x1", ylab="x2")
points(mu[,1], mu[,2], col="darkgreen")
persp(xx1, xx2, yyA, col="lightblue", xlab="x1", ylab="x2", zlab="f(x1, x2)")
contour(xx1, xx2, yyB, col="red", xlab="x1", ylab="x2")
points(mu[,1], mu[,2], col="darkgreen")
persp(xx1, xx2, yyB, col="lightblue", xlab="x1", ylab="x2", zlab="f(x1, x2)", phi=30, theta=30)

par(mfrow=c(2, 2), bty="n")
contour(xx1, xx2, yyC, col="red", xlab="x1", ylab="x2")
points(mu[,1], mu[,2], col="darkgreen")
persp(xx1, xx2, yyC, col="lightblue", xlab="x1", ylab="x2", zlab="f(x1, x2)")
contour(xx1, xx2, yyD, col="red", xlab="x1", ylab="x2")
points(mu[,1], mu[,2], col="darkgreen")
persp(xx1, xx2, yyD, col="lightblue", xlab="x1", ylab="x2", zlab="f(x1, x2)", phi=30, theta=30)
}
\keyword{distribution}
\keyword{multivariate}
