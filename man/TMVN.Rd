\name{TMVN}
\alias{TMVN}
\alias{rTMVN}
\title{
  Truncated multivariate normal distribution
}
\description{
  Random generation for the truncated multivariate normal distribution.
  The mean and covariance matrix of the original multivariate normal distribution
  are \code{mean} and \code{Sigma}. Truncation limits are given by
  \code{a}, \code{b}, type of truncation is given by \code{trunc}.

  This function uses a Gibbs algorithm to produce a Markov chain whose
  stationary distribution is the targeted truncated multivariate normal
  distribution, see Geweke (1991) for more details. Be aware that the
  sampled values are not i.i.d.!
}
\usage{
rTMVN(n, mean=c(0, 0), Sigma=diag(2), a, b, trunc, xinit)
}
\arguments{
  \item{mean}{a numeric vector of the mean of the original multivariate
    normal distribution.}
  \item{Sigma}{covariance matrix of the original multivariate normal
    distribution.}
  \item{a}{a numeric vector of the same length as \code{mean} of
    truncation limits 1.}
  \item{b}{a numeric vector of the same length as \code{mean} of
    truncation limits 2.}
  \item{trunc}{a numeric vector of the same length as \code{mean}
    describing the type of truncation in each margin.
    \describe{
      \item{\code{trunc}=0}{normal distribution is truncated on the
	interval \eqn{(a,\,\infty).}{(a, Infty)}. Value of \eqn{b} is
	ignored.}
      \item{\code{trunc}=1}{degenerated normal distribution, all values
	are with probability 1 equal to \eqn{a}, \eqn{b} is ignored.}
      \item{\code{trunc}=2}{normal distribution is truncated on the
	interval \eqn{(-\infty,\,a).}{(-Infty, a).} Value of \eqn{b} is ignored.}
      \item{\code{trunc}=3}{normal distribution is truncated on the
        interval \eqn{(a,\,b).}{(a,\,b).}}
      \item{\code{trunc}=4}{there is no truncation, values of \eqn{a}
	and \eqn{b} are ignored.}      
    }

    If \code{trunc} is not given, it is assumed that it is equal to
    4. Note that \code{a}, \code{b} and \code{trunc}
    must have the same length, with exception that 
    \code{b} does not have to be supplied if all \code{trunc} values 0,
    1, 2 or 4.
  }
  \item{xinit}{a numeric vector of the same length as \code{mean} with
    the initial value for the Gibbs sampler. If it is not supplied, the
    function determines itself the initial value.}
  \item{n}{number of observations to be sampled.}
}
\value{
  A matrix with the sampled values (Markov chain) in rows.
}
\references{
  Geweke, J. (1991).
  Efficient simulation from the multivariate normal and
  Student-t distributions subject to linear constraints 
  and the evaluation of constraint probabilities.
  \emph{Computer Sciences and Statistics}, \bold{23}, 571--578.
}
\seealso{
  \code{\link{rTNorm}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

exam2 <- function(n, mu, sigma, rho, a, b, trunc)
{
  Sigma <- matrix(c(sigma[1]^2, rho*sigma[1]*sigma[2], rho*sigma[1]*sigma[2], sigma[2]^2), nrow=2)
  x <- rTMVN(n, mean=mu, Sigma=Sigma, a=a, b=b, trunc=trunc)
  x1.gr <- seq(mu[1]-3.5*sigma[1], mu[1]+3.5*sigma[1], length=100)
  x2.gr <- seq(mu[2]-3.5*sigma[2], mu[2]+3.5*sigma[2], length=100)    
  z <- cbind(rep(x1.gr, 100), rep(x2.gr, each=100))
  dens.z <- matrix(dMVN(z, mean=mu, Sigma=Sigma), ncol=100)

  MEAN <- round(apply(x, 2, mean), 3)
  SIGMA <- var(x)
  SD <- sqrt(diag(SIGMA))
  RHO <- round(SIGMA[1,2]/(SD[1]*SD[2]), 3)
  SD <- round(SD, 3)
  
  layout(matrix(c(0,1,1,0, 2,2,3,3), nrow=2, byrow=TRUE))
  contour(x1.gr, x2.gr, dens.z, col="darkblue", xlab="x[1]", ylab="x[2]")
  points(x[,1], x[,2], col="red")
  title(sub=paste("Sample mean = ", MEAN[1], ", ", MEAN[2], ",  Sample SD = ", SD[1], ", ", SD[2], ",  Sample rho = ", RHO, sep=""))
  plot(1:n, x[,1], type="l", xlab="Iteration", ylab="x[1]", col="darkgreen")
  plot(1:n, x[,2], type="l", xlab="Iteration", ylab="x[2]", col="darkgreen")

  return(x)  
}  

x1 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0, a=c(-6, -9), b=c(4, 11), trunc=c(3, 3))
x2 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-6, -9), b=c(4, 11), trunc=c(3, 3))
x3 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-100, -100), b=c(100, 100), trunc=c(3, 3))
x4 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=-0.7, a=c(-6, -9), b=c(4, 11), trunc=c(3, 3))
x5 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=-0.9, a=c(-6, -9), b=c(4, 11), trunc=c(3, 3))

x6 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(0, 0), trunc=c(0, 2))
x7 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-1, 1), trunc=c(0, 2))
x8 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-1, 1), trunc=c(1, 2))
x9 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-1.5, 0.5), b=c(-0.5, 1.5), trunc=c(3, 3))
x10 <- exam2(1000, mu=c(-1, 1), sigma=c(1, sqrt(2)), rho=0.7, a=c(-1.5, 0.5), b=c(-0.5, 1.5), trunc=c(4, 3))
}
\keyword{distribution}
\keyword{multivariate}
