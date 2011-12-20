\name{GLMM_MCMCinit.alpha}
\alias{GLMM_MCMCinit.alpha}
\title{
  Handle init.alpha or init2.alpha argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCinit.alpha(init.alpha, lalpha, ialpha, number="")
}
\arguments{
  \item{init.alpha}{see argument \code{init.alpha} or \code{init2.alpha} of
    \code{\link{GLMM_MCMC}} function}
  \item{lalpha}{see output element \code{lalpha} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{ialpha}{a numeric vector of initial values of fixed effects.
    See output element \code{ialpha} of \code{\link{GLMM_MCMCifit}}
    function}
  \item{number}{label of the initial values (used only in error messages)}  
}
\value{
  A list with the same structure as complete \code{init.alpha} argument
  of the \code{\link{GLMM_MCMC}} function
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
