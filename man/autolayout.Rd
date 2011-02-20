\name{autolayout}
\alias{autolayout}
\title{
  Automatic layout for several plots in one figure
}
\description{
  Returns a matrix which can be used in \code{\link[graphics]{layout}}
  function as its \code{mat} argument.
}
\usage{
autolayout(np)
}
\arguments{
\item{np}{number of plots that are to be produced on 1 figure}
}
\value{
  A matrix.
}
\seealso{
  \code{\link{par}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
autolayout(10)
}
\keyword{dplot}
