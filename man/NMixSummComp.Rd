\name{NMixSummComp}
\alias{NMixSummComp}
\title{
  Summary for the mixture components
}
\description{
  This function returns basic posterior summary for (re-labeled) mixture
  components in a model with fixed number of components fitted with
  \code{\link{NMixMCMC}} function. The summary also takes into account
  possible scaling and shifting of the data (see argument \code{scale}
  in \code{\link{NMixMCMC}} function).

  Note that even though the mixture components are re-labeled before the
  summary is computed to achieve identifiability, posterior summaries of
  individual mixture means and variances are not always the quantity we
  would like to see. For density estimation, usually posterior
  predictive density (\code{\link{NMixPredDensMarg}},
  \code{\link{NMixPredDensJoint2}}) is usually the right stuff one
  should be interested in.
}
\usage{
NMixSummComp(x)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC}}
}
\value{
  Invisible \code{x}. The rest is printed on output device.
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
