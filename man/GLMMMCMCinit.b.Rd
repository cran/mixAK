\name{GLMM_MCMCinit.b}
\alias{GLMM_MCMCinit.b}
\title{
  Handle init.b or init2.b argument of GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  
  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCinit.b(init.b, prior.b, scale.b,
                id, dimb, LTb, naamLTb,
                I, ibMat, iEranefVec, iSDranefVec, number="")
}
\arguments{
  \item{init.b}{see argument \code{init.b} or \code{init2.b} of
    \code{\link{GLMM_MCMC}} function}
  \item{prior.b}{see argument \code{prior.b} of \code{\link{GLMM_MCMC}}
    function}
  \item{scale.b}{see argument \code{scale.b} of \code{\link{GLMM_MCMC}}
    function}
  \item{id}{}
  \item{dimb}{dimension of random effects (random intercepts included)}
  \item{LTb}{length of lower triangle of matrix dimb x dimb}
  \item{naamLTb}{names (indices) for elements of a lower triangle of matrix dimb x dimb}
  \item{I}{number of subjects (grouped observations) in the original data (before removing NA's)}
  \item{ibMat}{matrix with initial values of random effects}
  \item{iEranefVec}{vector with initial values of overall means of
    random effects}
  \item{iSDranefVec}{vector with initial values of overall standard
    deviations of random effects}  
  \item{number}{label of the initial values (used only in error messages)}  
}
\value{
  A list with the same structure as complete \code{init.b} argument
  of the \code{\link{GLMM_MCMC}} function
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{internal}
