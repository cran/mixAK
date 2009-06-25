\name{NMixChainsDerived}
\alias{NMixChainsDerived}
\title{
  Create MCMC chains derived from previously sampled values
}
\description{
  Currently, this function creates chains for marginal means
  of exp(data) from previously sampled values (see \code{\link{NMixMCMC}}).
  This is useful in survival context when a density
  of \eqn{Y=\log(T)}{Y=log(T)} is modelled using the function
  \code{\link{NMixMCMC}} and we are interested in inference
  on \eqn{\mbox{E}T = \mbox{E}\exp(Y)}{ET = Eexp(Y)}.
}
\usage{
NMixChainsDerived(object)
}
\arguments{
  \item{object}{an object of class \code{NMixMCMC} or \code{NMixMCMClist}}
}
\value{
  An object of the same class as argument \code{object}. When
  \code{object} was of class \code{NMixMCMC}, the resulting object
  contains additionally the following components:
  \item{chains.derived}{a \code{data.frame} with columns labeled
    \code{expy.Mean.1}, \ldots, \code{expy.Mean.p} containing the
    sampled values of \eqn{\mbox{E}\exp(Y_1)}{Eexp(Y[1])}, \ldots,
    \eqn{\mbox{E}\exp(Y_p)}{Eexp(Y[p])}.}
  \item{summ.expy.Mean}{posterior summary statistics for \eqn{\mbox{E}\exp(Y_1)}{Eexp(Y[1])}, \ldots,
    \eqn{\mbox{E}\exp(Y_p)}{Eexp(Y[p])}.}
   
  When \code{object} was of the class \code{NMixMCMClist} then each of
  its components (chains) is augmented by new components
  \code{chains.derived} and \code{summ.expy.Mean}.   
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
