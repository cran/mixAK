\name{Galaxy}
\alias{Galaxy}
\docType{data}
\title{Velocities of distant galaxies}
\description{
  Velocities (in km/sec) of 82 distant galaxies, diverging from our own galaxy.

  The dataset was first described by Roeder (1990) and subsequently
  analysed under different mixture models by several researchers
  including Escobar and West (1995) and Phillips and Smith (1996).
  Richardson and Green (1997) used a normal mixture estimated using
  reversible jump MCMC to estimate the distribution of the velocities.

  \bold{REMARK:} 78th observation is here 26.96 whereas it should be
  26.69 (see \code{\link[MASS]{galaxies}}). A value of 26.96 used in
  Richardson and Green (1997) is kept here.
  
}
\usage{data(Galaxy)}
\format{
  A numeric vector with observed values.
}
\source{
  \url{http://www.stats.bris.ac.uk/~peter/mixdata}
}
\references{
Escobar, M. D. and West, M. (1995).
Bayesian density estimation and inference using mixtures.
\emph{Journal of the American Statistical Association}, \bold{90}, 577--588.

Phillips, D. B. and Smith, A. F. M. (1996).
Bayesian model comparison via jump diffusions.
In \emph{Practical Markov Chain Monte Carlo},
eds: W. R. Gilks, S. Richardson, and D. J. Spiegelhalter, ch. 13, pp. 215-239.
London: Chapman and Hall.

Richardson, S. and Green, P. J. (1997).
On Bayesian analysis of mixtures with unknown number of components (with Discussion).
\emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.

Roeder, K. (1990).
Density estimation with confidence sets exemplified by superclusters and voids in the galaxies.
\emph{Journal of the American Statistical Association}, \bold{85}, 617--624.
}
\seealso{\code{\link[MASS]{galaxies}}.}
\examples{
data(Galaxy)
summary(Galaxy)
}
\keyword{datasets}
