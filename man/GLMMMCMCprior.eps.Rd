\name{GLMM_MCMCprior.eps}
\alias{GLMM_MCMCprior.eps}
\title{
  Handle prior.eps argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCprior.eps(prior.eps, Rc, isigma, is.sigma)
}
\arguments{
  \item{prior.eps}{see argument\code{prior.eps} of
    \code{\link{GLMM_MCMC}} function}
  \item{Rc}{see output element \code{Rc} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{isigma}{see output element \code{isigma} of \code{\link{GLMM_MCMCifit}}
    function. It is used to determine a default value of \code{prior.eps$h}.
    Defaults to a vector of ones if missing}
  \item{is.sigma}{see output element \code{is.sigma} of \code{\link{GLMM_MCMCifit}}
    function}  
}
\value{
  A list with the following components:
  \item{prior.eps}{value of \code{prior.eps} with added default values
    where necessary}
  \item{CpriorDouble_eps}{a numeric vector to be passed to C++}
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
