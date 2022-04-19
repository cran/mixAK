\name{summaryDiff}
\alias{summaryDiff}
\title{
  Posterior summary statistics for a difference of two quantities
}
\description{
  It calculates (posterior) summary statistics for a difference of two quantities supplied as (MCMC) samples
  Within \code{mixAK} package it is primarily used to calculate posterior summary for the difference of the
  deviances of two competing models.
}
\usage{
summaryDiff(x, y, prob=c(0.025, 0.5, 0.975), cut=c(-2*log(9), 0), na.rm=TRUE)
}
\arguments{
  \item{x}{a numeric vector with the sample of the first quantity.}
  \item{y}{a numeric vector with the sample of the second quantity to be subtracted from \code{x}.}
  \item{prob}{a numeric vector of probabilities for quantiles to be calculated from the sample of differences.}
  \item{cut}{numeric value(s) which specify the cutoff(s) we are interested in estimating
    \eqn{\mbox{P}(x - y < \mbox{cut})}{P(x - y < cut)} from the sample. The default values are motivated by
    the arguments given in Section 4 of Aitkin, Liu and Chadwick (2009)
    and in Section 7.5 of Aitkin (2010).
  }
  \item{na.rm}{logical indicating on how to handle \code{NA}'s.}
}
\value{
  A list with the components
  \item{summary}{a named vector with the (posterior) summary statistics
    based on the differences.}
  \item{Pcut}{estimated (posterior) probabilities that the difference
    lies below the \code{cut} values.}  
}
\references{
  Aitkin, M. (2010).
  \emph{Statistical Inference: An Integrated Bayesian/Likelihood Approach.}
  Boca Raton: CRC Press.
  
  Aitkin, M., Liu, C. C., and Chadwick, T. (2009).
  Bayesian model comparison and model averaging for small-area
  estimation.
  \emph{Annals of Applied Statistics}, \bold{3}, 199-221.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
set.seed(16336886)
x <- runif(100, 0, 100)
y <- runif(100, 0, 100)
sdiff <- summaryDiff(x, y)
}
\keyword{univar}
\keyword{htest}
