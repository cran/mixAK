\name{GLMM_MCMCinit.eps}
\alias{GLMM_MCMCinit.eps}
\title{
  Handle init.eps or init2.eps argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCinit.eps(init.eps, prior.eps, Rc, isigma, number="")
}
\arguments{
  \item{init.eps}{see argument \code{init.eps} or \code{init2.eps} of \code{\link{GLMM_MCMC}} function}
  \item{prior.eps}{see argument \code{prior.eps} of \code{\link{GLMM_MCMC}} function}  
  \item{Rc}{see output element \code{Rc} of \code{\link{GLMM_MCMCdata}} function}
  \item{isigma}{a numeric vector of length \code{Rc} with initial values of residual standard deviations.}
  \item{number}{label of the initial values (used only in error messages)}  
}
\value{
  A list with the same structure as complete \code{init.eps} argument
  of the \code{\link{GLMM_MCMC}} function
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
