\name{NMixMCMCinitr}
\alias{NMixMCMCinitr}
\title{
  Initial component allocations for the NMixMCMC function 
}
\description{
  This is a help function for \code{\link{NMixMCMC}} function. If
  \code{initr} is not given, it calculates reasonable initial values for
  component allocations. If \code{initr} argument is given then it is
  checked for consistency and formatted on output.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
NMixMCMCinitr(z, K, w, mu, Sigma, p, n, initr)
}
\arguments{
  \item{z}{a vector (if \eqn{p=1}) or a \eqn{n\times p}{n x p} matrix
    with initial values of scaled observations}
  \item{K}{a number of mixture components}
  \item{w}{a vector with (initial) mixture weights}
  \item{mu}{a vector of (initial) mixture means (if \eqn{p=1}) or a
    \eqn{K\times p}{K x p} matrix with (initial) mixture means in rows
    (if \eqn{p>1})}
  \item{Sigma}{a vector of (initial) mixture variances (if \eqn{p=1}) or
    a \eqn{K\,p\times p}{K*p x p} matrix with (initial) mixture
    covariance matrices put below each other}
  \item{p}{dimension of the data}  
  \item{n}{number of observations}
  \item{initr}{a vector of initial component allocations to be checked
    for consistency. If not given then reasonable initials are generated.}
}
\value{
A vector with reasonable initial component allocations.
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
