\name{NMixMCMCinity}
\alias{NMixMCMCinity}
\title{
  Initial values of censored observations for the NMixMCMC function 
}
\description{
  This is a help function for \code{\link{NMixMCMC}} function. If
  \code{inity} is not given, it calculates reasonable initial values for
  censored observations. If \code{inity} argument is given then it is
  checked for consistency and formatted on output.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
NMixMCMCinity(y0, y1, censor, sd.init,
   are.Censored, are.Right, are.Exact, are.Left, are.Interval,
   p, n, inity, random=FALSE)
}
\arguments{
  \item{y0}{see output from \code{\link{NMixMCMCdata}}}
  \item{y1}{see output from \code{\link{NMixMCMCdata}}}
  \item{censor}{see output from \code{\link{NMixMCMCdata}}}  
  \item{sd.init}{a vector of length \eqn{p} with initial values for
    overall standard deviations in each margin}
  \item{are.Censored}{see output from \code{\link{NMixMCMCdata}}}
  \item{are.Right}{see output from \code{\link{NMixMCMCdata}}} 
  \item{are.Exact}{see output from \code{\link{NMixMCMCdata}}}
  \item{are.Left}{see output from \code{\link{NMixMCMCdata}}}
  \item{are.Interval}{see output from \code{\link{NMixMCMCdata}}}      
  \item{p}{dimension of the data}  
  \item{n}{number of observations}
  \item{inity}{a vector (if \eqn{p=1}) or a \eqn{n\times p}{n x p} matrix
      (if \eqn{p \geq 1}{p >= 1}) of initial values of censored observations to be checked
      for consistency. If not given then reasonable initials are
      generated.}
  \item{random}{logical value. If \code{TRUE} then some randomness is
    used when generating initial values.}  
}
\value{
A \eqn{n\times p}{n x p} matrix with reasonable initial values for
censored observations.
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
