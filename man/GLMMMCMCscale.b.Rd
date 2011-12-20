\name{GLMM_MCMCscale.b}
\alias{GLMM_MCMCscale.b}
\title{
  Handle scale.b argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCscale.b(scale.b, dimb, iEranefVec, iSDranefVec)
}
\arguments{
  \item{scale.b}{see argument\code{prior.eps} of
    \code{\link{GLMM_MCMC}} function}
  \item{dimb}{see output element \code{dimb} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{iEranefVec}{see output element \code{iEranefVec} of \code{\link{GLMM_MCMCifit}}
    function. It is used to determine a default value of \code{scale.b$shift}.
    Defaults to a vector of zeros if missing}
  \item{iSDranefVec}{see output element \code{iSDranefVec} of \code{\link{GLMM_MCMCifit}}
    function. It is used to determine a default value of \code{scale.b$scale}.
    Defaults to a vector of ones if missing}  
}
\value{
  Value of \code{scale.b} with added default values where necessary.
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}

