\name{Enzyme}
\alias{Enzyme}
\docType{data}
\title{Enzymatic activity in the blood}
\description{
  Enzymatic activity in the blood, for an enzyme involved in the
  metabolism of carcinogenic substances, among a group of 245 unrelated individuals.

  Bechtel et al. (1993) identified a mixture of two skewed distributions
  by using maximum-likelihood estimation.
  Richardson and Green (1997) used a normal mixture estimated using
  reversible jump MCMC to estimate the distribution of the enzymatic activity.
   
}
\usage{data(Enzyme)}
\format{
  A numeric vector with observed values.
}
\source{
  Originally from http://www.stats.bris.ac.uk/~peter/mixdata/
}
\references{
\eqn{\mbox{Bechtel, Y. C., Bona\"{\i}ti-Pelli\'e, C., Poisson, N., Magnette, J., and Bechtel, P. R.}}{%
     Bechtel, Y. C., Bonaiti-Pellie, C., Poisson, N., Magnette, J., and Bechtel, P. R.}(1993).  
A population and family study of N-acetyltransferase using caffeine urinary metabolites.
\emph{Clinical Pharmacology and Therapeutics}, \bold{54}, 134--141.
  
Richardson, S. and Green, P. J. (1997).
On Bayesian analysis of mixtures with unknown number of components (with Discussion).
\emph{Journal of the Royal Statistical Society, Series B}, \bold{59}, 731--792.
}
\examples{
data(Enzyme)
summary(Enzyme)
}
\keyword{datasets}
