\name{NMixMCMCdata}
\alias{NMixMCMCdata}
\title{
  Data manipulation for the NMixMCMC function 
}
\description{
  This is a help function for \code{\link{NMixMCMC}} function.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
NMixMCMCdata(y0, y1, censor)
}
\arguments{
  \item{y0}{see \code{\link{NMixMCMC}} function.}
  \item{y1}{see \code{\link{NMixMCMC}} function.}
  \item{censor}{see \code{\link{NMixMCMC}} function.}  
}
\value{
  A list with the following components:
  \item{y0}{worked out value of \code{y0}}
  \item{y1}{worked out value of \code{y1}}
  \item{censor}{worked out value of \code{censor}}
  \item{n}{number of observations}
  \item{p}{dimension of the response}
  \item{LTp}{length of the lower triangle of the \eqn{p\times p}{p x p} matrix}
  \item{naamLTp}{a character vector of length \code{LTp} with names
    (indices) for elements of a lower triangle of the \eqn{p\times p}{p x p} matrix}
  \item{are.Censored}{logical value which indicates whether there are
    censored observations in the data}
  \item{are.Right}{logical value which indicates whether there are
    right-censored observations in the data}
  \item{are.Exact}{logical value which indicates whether there are
    exact observations in the data}
  \item{are.Left}{logical value which indicates whether there are
    left-censored observations in the data}  
  \item{are.Interval}{logical value which indicates whether there are
    interval-censored observations in the data}
  \item{is.Interval}{a logical vector of length \eqn{n} which indicates
    whether a specific observation is interval-censored}
}
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
