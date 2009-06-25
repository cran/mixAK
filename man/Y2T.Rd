\name{Y2T}
\alias{Y2T}
\alias{Y2T.NMixPredDensMarg}
\alias{Y2T.NMixPredCDFMarg}
\title{
  Transform fitted distribution of Y=log(T) into distribution of T
}
\description{
  This method transforms fitted distributionof \eqn{Y=\log(T)}{Y=log(T)}
  into distribution of \eqn{T}. 
}
\usage{
Y2T(x, \dots)

\method{Y2T}{NMixPredDensMarg}(x, \dots)

\method{Y2T}{NMixPredCDFMarg}(x, \dots)
}
\arguments{
  \item{x}{an object of appropriate class.}
  \item{\dots}{optional additional arguments.}  
}
\value{
  An object of the same class as argument \code{x}.
}
\seealso{
  \code{\link{NMixPredDensMarg}}, \code{\link{NMixPredCDFMarg}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{smooth}
