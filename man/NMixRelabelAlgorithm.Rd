\name{NMixRelabelAlgorithm}
\alias{NMixRelabelAlgorithm}
\title{
  Argument manipulation for the NMixRelabel functions 
}
\description{
  This is a help function for \code{\link{NMixRelabel}} functions.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
NMixRelabelAlgorithm(type=c("mean", "weight", "stephens"), par, dim)
}
\arguments{
  \item{type}{see \code{\link{NMixRelabel}} function.}
  \item{par}{see \code{\link{NMixRelabel}} function.}  
  \item{dim}{dimension of the data}  
}
\value{
  A list with the following components:
  \item{relabel}{a list with the same structure as \code{relabel}
    component of objects of class \code{NMixMCMC} or \code{GLMM_MCMC}.}
  \item{Ctype}{integer to be passed to C++ re-labelling functions.}
  \item{iparam}{integer vector to be passed to C++ re-labelling functions.}  
}  
\seealso{
  \code{\link{NMixRelabel}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{internal}
