\name{NMixPlugDA}
\alias{NMixPlugDA}
\title{
  Discriminant analysis based on plug-in estimates from the mixture model
}
\description{
  It performs discriminant analysis based on 
  posterior summary for (re-labeled) mixture
  components in a model with fixed number of components fitted with
  \code{\link{NMixMCMC}} function.
}
\usage{
NMixPlugDA(object, y)
}
\arguments{
  \item{object}{an object of class \code{NMixMCMC}}
  \item{y}{vector, matrix or data frame with observations to be clustered}
}
\value{
  A \code{data.frame} with 
  columns labeled \code{prob1},\dots, \code{probp} giving
  plug-in estimates of probabilities of belonging to each component and a column
  labeled \code{component} giving the index of the component with the
  highest component probability.  
}
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{NMixPredDA}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
\keyword{cluster}
