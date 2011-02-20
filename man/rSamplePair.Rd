\name{rSamplePair}
\alias{rSamplePair}
\title{
  Sample a pair (with replacement)
}
\description{
  For given \eqn{K}, the function samples with replacement from a
  uniform distribution on a set of pairs
  \eqn{(1,2), (1,3), \ldots, (1,K), (2,3), \ldots, (2,K), \ldots, (K-1,K).}
}
\usage{
  rSamplePair(n, K)
}
\arguments{
  \item{n}{number of pairs to sample.}
  \item{K}{a numeric value which determines \eqn{K} (see above).}
}
\value{
  A two-component numeric vector for \code{n}\eqn{=2} or a matrix with 2
  columns with sampled pairs in rows for \code{n}\eqn{>2.}
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
rSamplePair(n=1, K=2)
rSamplePair(n=10, K=2)

rSamplePair(n=1, K=3)
rSamplePair(n=10, K=3)

rSamplePair(n=1, K=4)
rSamplePair(n=10, K=4)
}
\keyword{distribution}
