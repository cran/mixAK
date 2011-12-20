\name{NMixMCMCwrapper}
\alias{NMixMCMCwrapper}
\title{
  Wrapper to the NMixMCMC main simulation.
}
\description{
  This is wrapper to the NMixMCMC main simulation which allows
  vectorized evaluation and possibly parallel computation.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
NMixMCMCwrapper(chain=1,
                scale, prior, inits, Cpar, RJMCMC, CRJMCMC,
                actionAll, nMCMC, keep.chains, PED,
                dens.zero)
}
\arguments{
  \item{chain}{identification of the chain sampled in a particular call
    of this function, usually number like 1, 2, \ldots}
  \item{Cpar}{a list with the following components
    \describe{
      \item{z0}{\eqn{n\times p}{n x p} matrix with shifted and scaled main limits
        of observed intervals.}
      \item{z0}{\eqn{n\times p}{n x p} matrix with shifted and scaled upper
        limits of observed intervals.}
      \item{censor}{\eqn{n\times p}{n x p} matrix with censoring indicators.}
      \item{p}{dimension of the response.}
      \item{n}{number of observations.}
      \item{Cinteger}{a numeric vector with integer prior parameters.}
      \item{Cdouble}{a numeric vector with double precission prior parameters.}      
    }  
  }  
  \item{scale}{a list specifying how to scale the data before running
    MCMC. See argument \code{scale} in \code{\link{NMixMCMC}}
  }
  \item{prior}{a list specifying prior hyperparameters. See argument
    \code{prior} in \code{\link{NMixMCMC}}.
  }
  \item{inits}{a list of length at least \code{chain}. Its
    \code{chain}-th component is used. Each component of the list should
    have the structure of \code{init} argument of function
    \code{\link{NMixMCMC}}.
  }
  \item{RJMCMC}{a list specifying parameters for RJ-MCMC.
    See argument \code{RJMCMC} in \code{\link{NMixMCMC}}}
  \item{CRJMCMC}{a numeric vector with parameters for RJ-MCMC.}
  \item{actionAll}{argument for underlying C++ function.}
  \item{nMCMC}{vector giving the length of MCMC etc.}
  \item{keep.chains}{logical. If \code{FALSE}, only summary statistics
    are returned in the resulting object. This might be useful in the
    model searching step to save some memory.}
  \item{PED}{a logical value which indicates whether the penalized
    expected deviance (see Plummer, 2008 for more details)
    will be computed (which requires two parallel
    chains). Even if \code{keep.chains} is \code{FALSE}, it is necessary
    to keep (for a while) at least some chains to compute PED.
  }
  \item{dens.zero}{small number (1e-300) to determine whether the
    contribution to the deviance (\eqn{-\log}{-log} density) is equal to
    infinity. Such values are trimmed when computing expected deviance.}
}  
\value{
  A list having almost the same components as object returned by
  \code{\link{NMixMCMC}} function.
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
