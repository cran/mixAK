\name{NMixClust}
\alias{NMixClust}
\title{
  Clustering based on a mixture model
}
\description{
  It performs clustering based on posterior summary for (re-labeled) mixture
  components in a model with fixed number of components fitted with
  \code{\link{NMixMCMC}} function.
}
\usage{
NMixClust(object, y)
}
\arguments{
  \item{object}{an object of class \code{NMixMCMC}}
  \item{y}{vector, matrix or data frame with observations to be clustered}
}
\value{
  A \code{data.frame} with observations in the first columns and
  additionally columns labeled \code{prob1},\dots, \code{probp} giving
  posterior probabilities of belonging to each component and a column
  labeled \code{component} giving the index of the component with the
  highest posterior probability.  
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
