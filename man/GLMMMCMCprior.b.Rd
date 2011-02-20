\name{GLMM_MCMCprior.b}
\alias{GLMM_MCMCprior.b}
\title{
  Handle prior.eps argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCprior.b(prior.b, scale.b, dimb, iEranefVec, iSDranefVec)
}
\arguments{
  \item{prior.b}{see argument\code{prior.b} of
    \code{\link{GLMM_MCMC}} function}
  \item{scale.b}{see argument\code{scale.b} of
    \code{\link{GLMM_MCMC}} function}
  \item{dimb}{see output element \code{dimb} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{iEranefVec}{see output element \code{iEranefVec} of \code{\link{GLMM_MCMCifit}}
    function. Defaults to a vector of zeros if missing}
  \item{iSDranefVec}{see output element \code{iSDranefVec} of \code{\link{GLMM_MCMCifit}}
    function. Defaults to a vector of ones if missing}  
}
\value{
  A list with the following components:
  \item{prior.b}{value of \code{prior.b} with added default values
    where necessary}
  \item{CpriorInt_b}{a numeric vector to be passed to C++}  
  \item{CpriorDouble_b}{a numeric vector to be passed to C++}
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
