\name{Acidity}
\alias{Acidity}
\docType{data}
\title{Acidity index of lakes in North-Central Wisconsin}
\description{
  Acidity index measured in a sample of 155 lakes in North-Central
  Wisconsin.

  Crawford et al. (1992) and Crawford (1994) analyzed these data as a
  mixture of normal distributions on the log-scale. Richardson and Green
  (1997) used a normal mixture estimated using reversible jump MCMC.
}
\usage{data(Acidity)}
\format{
  A numeric vector with observed values.
}
\source{
  \url{http://www.stats.bris.ac.uk/~peter/mixdata}
}
\references{
Crawford, S. L. (1994).
An application of the Laplace method to finite mixture distribution.
\emph{Journal of the American Statistical Association}, \bold{89}, 259--267.

Crawford, S. L., DeGroot, M. H., Kadane, J. B., and Small, M. J. (1994).
Modeling lake chemistry distributions: Approximate Bayesian methods for estimating a finite mixture model.
\emph{Technometrics}, \bold{34}, 441--453.

Richardson, S. and Green, P. J. (1997).
On Bayesian analysis of mixtures with unknown number of components (with Discussion).
\emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.
}
\examples{
data(Acidity)
summary(Acidity)
}
\keyword{datasets}
