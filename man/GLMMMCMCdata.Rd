\name{GLMM_MCMCdata}
\alias{GLMM_MCMCdata}
\title{
  Data manipulation for the GLMM_MCMC function 
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCdata(y, dist, id, time, x, z, random.intercept)
}
\arguments{
  \item{y}{see \code{\link{GLMM_MCMC}} function.}
  \item{dist}{see \code{\link{GLMM_MCMC}} function.}
  \item{id}{see \code{\link{GLMM_MCMC}} function.}
  \item{time}{see \code{\link{GLMM_longitDA}} function.}  
  \item{x}{see \code{\link{GLMM_MCMC}} function.}
  \item{z}{see \code{\link{GLMM_MCMC}} function.}
  \item{random.intercept}{see \code{\link{GLMM_MCMC}} function.}  
}
\value{
  A list with the following components:
  \item{y}{worked out value of \code{y}}
  \item{dist}{worked out value of \code{dist}}  
  \item{id}{worked out value of \code{id}}
  \item{time}{worked out value of \code{time}}  
  \item{x}{worked out value of \code{x}}
  \item{z}{worked out value of \code{z}}
  \item{random.intercept}{worked out value of \code{random.intercept}}  
  \item{R}{number of response variables}
  \item{Rc}{number of continuous response variables}
  \item{Rd}{number of discrete response variables}  
  \item{ndist}{numerical counterpart of 'dist' argument}
  \item{xempty}{logical vector of length R}
  \item{zempty}{logical vector of length R}
  \item{p}{numeric vector of length R}  
  \item{q}{numeric vector of length R}
  \item{fixed.intercept}{logical vector of length R}
  \item{CrandomIntcpt}{numerical counterpart of random.intercept}
  \item{CfixedIntcpt}{numerical counterpart of fixed.intercept}
  \item{dimb}{dimension of random effects (random intercepts included)}
  \item{LTb}{length of lower triangle of matrix dimb x dimb}
  \item{naamLTb}{names (indices) for elements of a lower triangle of matrix dimb x dimb}  
  \item{lalpha}{length of alpha vector (fixed intercepts included)}
  \item{p_fi}{numeric vector of length R}
  \item{q_ri}{numeric vector of length R}
  \item{name.response}{character vector}
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
