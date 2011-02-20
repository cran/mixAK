\name{SP2Rect}
\alias{SP2Rect}
\title{
  Creation of a rectangular symmetric matrix from its lower triangle
  (symmetric packed format)
}
\description{
  It creates a symmetric matrix from its lower triangle.
}
\usage{
  SP2Rect(LT, dim)
}
\arguments{
  \item{LT}{a numeric vector with the lower triangle (stored columnwise)
    of the matrix we want to reconstruct.}
  \item{dim}{number of rows and columns of a resulting matrix.}  
}
\value{
  A matrix.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
SP2Rect(3, dim=1)
SP2Rect(c(1, 0.5, 2), dim=2)
SP2Rect(c(1, 0.5, 0.25, 2, -0.5, 3), dim=3)
}
\keyword{array}
