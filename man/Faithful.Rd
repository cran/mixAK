\name{Faithful}
\alias{Faithful}
\docType{data}
\title{Old Faithful Geyser Data}
\description{
    Waiting time between eruptions and the duration of the eruption
    for the Old Faithful geyser in Yellowstone National Park, Wyoming,
    USA, version used in \eqn{\mbox{H\"ardle, W.}}{H&#228rdle, W.} (1991).
}
\details{
  There are many versions of this dataset around. Azzalini and Bowman
  (1990) use a more complete version.  
}  
\usage{data(Faithful)}
\format{
  A data frame with 272 observations on 2 variables.
  \describe{
    \item{\code{eruptions}}{eruption time in minutes}
    \item{\code{waiting}}{waiting time to the next eruption in minutes}    
  }  
}
\source{
  R package MASS
}
\references{
  \eqn{\mbox{H\"ardle, W.}}{Haerdle, W.} (1991).
  \emph{Smoothing Techniques with Implementation in S}.
  New York: Springer.

  Azzalini, A. and Bowman, A. W. (1990).
  A look at some data on the Old Faithful geyser.
  \emph{Applied Statistics}, \bold{39}, 357-365.
}
\seealso{\code{\link[MASS]{geyser}}.}
\examples{
data(Faithful)
summary(Faithful)
}
\keyword{datasets}
