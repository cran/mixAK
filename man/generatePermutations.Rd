\name{generatePermutations}
\alias{generatePermutations}
\title{
  Generate all permutations of (1, ..., K) 
}
\description{
  It generates a matrix containing all permutations of (1, ..., K).
}
\usage{
generatePermutations(K)
}
\arguments{
  \item{K}{integer value of \eqn{K}.}
}
\value{
  A matrix of dimension \eqn{K! \times K}{K! x K} with generated
  permutations in rows.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
generatePermutations(1)
generatePermutations(2)
generatePermutations(3)
generatePermutations(4)
}
\keyword{utilities}
\keyword{arith}
