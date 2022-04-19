\name{NMixSummComp}
\alias{NMixSummComp}
\alias{NMixSummComp.default}
\alias{NMixSummComp.NMixMCMC}
\alias{NMixSummComp.GLMM_MCMC}
\title{
  Summary for the mixture components
}
\description{
  This function returns basic posterior summary for (re-labeled) mixture
  components in a model with fixed number of components fitted with
  \code{\link{NMixMCMC}} or \code{\link{GLMM_MCMC}} function.
  The summary also takes into account
  possible scaling and shifting of the data (see argument \code{scale}
  in \code{\link{NMixMCMC}} function
  or argument \code{scale.b} in \code{\link{GLMM_MCMC}}).

  Note that even though the mixture components are re-labeled before the
  summary is computed to achieve some identifiability, posterior summaries of
  individual mixture means and variances are not always the quantity we
  would like to see. For density estimation, posterior
  predictive density (\code{\link{NMixPredDensMarg}},
  \code{\link{NMixPredDensJoint2}}) is usually the right stuff one
  should be interested in.
}
\usage{
NMixSummComp(x)

\method{NMixSummComp}{default}(x)

\method{NMixSummComp}{NMixMCMC}(x)

\method{NMixSummComp}{GLMM_MCMC}(x)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} or \code{GLMM_MCMC}}
}
\value{
  Invisible \code{x}. The rest is printed on output device.
}
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
