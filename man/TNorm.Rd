\name{TNorm}
\alias{TNorm}
\alias{rTNorm}
\alias{C_rTNorm1_R}
\title{
  Truncated normal distribution
}
\description{
  Random generation for the truncated normal distribution.
  The mean and standard deviation of the original normal distribution
  are \code{mean} and \code{sd}. Truncation limits are given by
  \code{a}, \code{b}, type of truncation is given by \code{trunc}.
}
\usage{
rTNorm(n, mean=0, sd=1, a, b, trunc)
}
\arguments{
  \item{mean}{mean (if common for all observations) or a vector of
    length \code{n} of means.}
  \item{sd}{standard deviation (if common for all observations) or a vector of
    length \code{n} of standard deviations.

    Note that \code{mean} and \code{sd}
    must have the same length, either 1 or \code{n}.}
  \item{a}{truncation limit 1 (if common for all observations) or a
    vector of length \code{n} of truncation limits 1.}
  \item{b}{truncation limit 2 (if common for all observations) or a
    vector of length \code{n} of truncation limits 2.}
  \item{trunc}{type of truncation (if common for all observations) or a
    vector of length \code{n} of types of truncation

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
    must have the same length, either 1 or \code{n} with exception that
    \code{b} does not have to be supplied if \code{trunc} is 0, 1, 2 or 4.    
  }    
  \item{n}{number of observations to be sampled.}
}
\value{
  A numeric vector with sampled values.
}
\references{
  Geweke, J. (1991).
  Efficient simulation from the multivariate normal and
  Student-t distributions subject to linear constraints 
  and the evaluation of constraint probabilities.
  \emph{Computer Sciences and Statistics}, \bold{23}, 571--578.
}
\seealso{
  \code{\link{rnorm}}, \code{\link{rTMVN}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

### Not truncated normal distribution
x1 <- rTNorm(1000, mean=10, sd=3)
c(mean(x1), sd(x1), range(x1))

### Truncation from left only
x2 <- rTNorm(1000, mean=10, sd=3, a=7, trunc=0)
c(mean(x2), sd(x2), range(x2))

### Degenerated normal distribution
x6 <- rTNorm(1000, mean=10, sd=3, a=13, trunc=1)
c(mean(x6), sd(x6), range(x6))

### Truncation from right only
x3 <- rTNorm(1000, mean=10, sd=3, a=13, trunc=2)
c(mean(x3), sd(x3), range(x3))

### Truncation from both sides
x4 <- rTNorm(1000, mean=10, sd=3, a=7, b=13, trunc=3)
c(mean(x4), sd(x4), range(x4))

x5 <- rTNorm(1000, mean=10, sd=3, a=5.5, b=14.5, trunc=3)
c(mean(x5), sd(x5), range(x5))

oldPar <- par(mfrow=c(2, 3))
hist(x1, main="N(10, 3^2)")
hist(x2, main="TN(10, 3^2, 7, Infty)")
hist(x6, main="TN(10, 3^2, 13, 13)")
hist(x3, main="TN(10, 3^2, -Infty, 13)")
hist(x4, main="TN(10, 3^2, 7, 13)")
hist(x5, main="TN(10, 3^2, 5.5, 14.5)")
par(oldPar)

### Different truncation limits
n <- 1000
a <- rnorm(n, -2, 1)
b <- a + rgamma(n, 1, 1)
trunc <- rep(c(0, 1, 2, 3, 4), each=n/5)
x7 <- rTNorm(n, mean=1, sd=2, a=a, b=b, trunc=trunc)
cbind(trunc, a, x7)[1:10,]
sum(x7[1:(n/5)] > a[1:(n/5)])      ## must be equal to n/5

cbind(trunc, a, x7)[201:210,]
sum(x7[(n/5+1):(2*n/5)] == a[(n/5+1):(2*n/5)])         ## must be equal to n/5

cbind(trunc, x7, a)[401:410,]
sum(x7[(2*n/5+1):(3*n/5)] < a[(2*n/5+1):(3*n/5)])      ## must be equal to n/5

cbind(trunc, a, x7, b)[601:610,]
sum(x7[(3*n/5+1):(4*n/5)] > a[(3*n/5+1):(4*n/5)])      ## must be equal to n/5
sum(x7[(3*n/5+1):(4*n/5)] < b[(3*n/5+1):(4*n/5)])      ## must be equal to n/5

cbind(trunc, x7)[801:810,]

### Different moments and truncation limits
n <- 1000
mu <- rnorm(n, 1, 0.2)
sigma <- 0.5 + rgamma(n, 1, 1)
a <- rnorm(n, -2, 1)
b <- a + rgamma(n, 1, 1)
trunc <- rep(c(0, 1, 2, 3, 4), each=n/5)
x8 <- rTNorm(n, mean=1, sd=2, a=a, b=b, trunc=trunc)

### Truncation from left only
### (extreme cases when we truncate to the area
###  where the original normal distribution has
###  almost zero probability)
x2b <- rTNorm(1000, mean=0, sd=1, a=7.9, trunc=0)
c(mean(x2b), sd(x2b), range(x2b))

x2c <- rTNorm(1000, mean=1, sd=2, a=16, trunc=0)
c(mean(x2c), sd(x2c), range(x2c))

### Truncation from right only (extreme cases)
x3b <- rTNorm(1000, mean=0, sd=1, a=-7.9, trunc=2)
c(mean(x3b), sd(x3b), range(x3b))

x3c <- rTNorm(1000, mean=1, sd=2, a=-13, trunc=2)
c(mean(x3c), sd(x3c), range(x3c))

### Truncation from both sides (extreme cases)
x4b <- rTNorm(1000, mean=0, sd=1, a=-9, b=-7.9, trunc=3)
c(mean(x4b), sd(x4b), range(x4b))

x4c <- rTNorm(1000, mean=0, sd=1, a=7.9, b=9, trunc=3)
c(mean(x4c), sd(x4c), range(x4c))
}
\keyword{distribution}

