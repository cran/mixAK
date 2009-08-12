\name{GLMM_MCMCprior.beta}
\alias{GLMM_MCMCprior.beta}
\title{
  Handle prior.beta argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCprior.beta(prior.beta, lbeta)
}
\arguments{
  \item{prior.beta}{see argument\code{prior.beta} of
    \code{\link{GLMM_MCMC}} function}
  \item{lbeta}{see output element \code{lbeta} of \code{\link{GLMM_MCMCdata}}
    function}  
}
\value{
  A list with the following components:
  \item{prior.beta}{value of \code{prior.beta} with added default values
    where necessary}
  \item{CpriorDouble_beta}{a numeric vector to be passed to C++}
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
