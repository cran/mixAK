\name{NMixPredDA}
\alias{NMixPredDA}
\alias{C_NMix_PredDA}
\title{
  Discriminant analysis based on MCMC output from the mixture model
}
\description{
  It performs discriminant analysis based on sampled (re-labeled) MCMC
  chains from the mixture model fitted with
  \code{\link{NMixMCMC}} function. Observations to be discriminated may
  be censored.

  Discrimination is based on posterior predictive probabilities of
  belonging to (re-labeled) mixture components.
}
\usage{
NMixPredDA(object, y0, y1, censor, inity, info)
}
\arguments{
  \item{object}{an object of class \code{NMixMCMC}}
  \item{y0}{vector, matrix or data frame with observations (or limits of
    censored-observations) to be clustered. See \code{\link{NMixMCMC}} for details.
    
    If \code{y0} is not given then the function discriminates original
    observations used to generate MCMC sample stored in \code{object}.
  }
  \item{y1}{vector, matrix or data frame with upper limits of
    interval-censored observations (if there are any). See
    \code{\link{NMixMCMC}} for details.
  }
  \item{censor}{vector, matrix or data frame with censoring indicators
    (if there are any censored observations). See
    \code{\link{NMixMCMC}} for details.
  }
  \item{inity}{optional vector, matrix or data frame with initial values
    of censored observations (if there are any censored observations)
  }
  \item{info}{number which specifies frequency used to re-display the
    iteration counter during the computation.
  }  
}
\value{
  A \code{data.frame} with 
  columns labeled \code{prob1},\dots, \code{probp} giving
  posterior predictive probabilities of belonging to each component and a column
  labeled \code{component} giving the index of the component with the
  highest component probability.  
}
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{NMixPlugDA}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
\keyword{cluster}
