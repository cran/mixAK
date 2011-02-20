\name{GLMM_MCMCprior.alpha}
\alias{GLMM_MCMCprior.alpha}
\title{
  Handle prior.alpha argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCprior.alpha(prior.alpha, lalpha)
}
\arguments{
  \item{prior.alpha}{see argument\code{prior.alpha} of
    \code{\link{GLMM_MCMC}} function}
  \item{lalpha}{see output element \code{lalpha} of \code{\link{GLMM_MCMCdata}}
    function}  
}
\value{
  A list with the following components:
  \item{prior.alpha}{value of \code{prior.alpha} with added default values
    where necessary}
  \item{CpriorDouble_alpha}{a numeric vector to be passed to C++}
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
