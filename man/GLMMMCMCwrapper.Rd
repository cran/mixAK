\name{GLMM_MCMCwrapper}
\alias{GLMM_MCMCwrapper}
\title{
  Wrapper to the GLMM_MCMC main simulation.
}
\description{
  This is wrapper to the GLMM_MCMC main simulation which allows
  vectorized evaluation and possibly parallel computation.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCwrapper(chain=1, data,
     prior.alpha, init.alpha,
     scale.b, prior.b, init.b,
     prior.eps, init.eps,
     Cpar, nMCMC, store, keep.chains)
}
\arguments{
  \item{chain}{identification of the chain sampled in a particular call
    of this function, usually number like 1, 2, \ldots}
  \item{data}{}
  \item{prior.alpha}{}
  \item{init.alpha}{a list of \code{init.alpha} objects for all chains}
  \item{scale.b}{}
  \item{prior.b}{}
  \item{init.b}{a list of \code{init.b} objects for all chains}
  \item{prior.eps}{}
  \item{init.eps}{a list of \code{init.eps} objects for all chains}
  \item{Cpar}{}
  \item{nMCMC}{}
  \item{store}{}
  \item{keep.chains}{}      
}  
\value{
  A list having almost the same components as object returned by
  \code{\link{GLMM_MCMC}} function.
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
