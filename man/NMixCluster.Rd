\name{NMixCluster}
\alias{NMixCluster}
\alias{NMixCluster.default}
\alias{NMixCluster.GLMM_MCMC}
\title{
  Clustering based on the MCMC output of the mixture model
}
\description{
  TO BE ADDED.
  
  This function only works for models with a fixed number of mixture components.
}
\usage{
NMixCluster(object, \dots)

\method{NMixCluster}{default}(object, \dots)

\method{NMixCluster}{GLMM_MCMC}(object,
   prob = c("poster.comp.prob", "quant.comp.prob", "poster.comp.prob_b",
            "quant.comp.prob_b", "poster.comp.prob_u"),
   pquant = 0.5, HPD = FALSE, pHPD = 0.95, pthresh = -1, unclass.na = FALSE, \dots)
}
\arguments{
  \item{object}{an object of apropriate class.}
  \item{prob}{character string which identifies estimates of the
    component probabilities to be used for clustering.}
  \item{pquant}{when \code{prob} is either \dQuote{quant.comp.prob} or
    \dQuote{quant.comp.prob_b}, argument \code{pquant} is the
    probability of the quantile of the component probabilities to be
    used for clustering.}
  \item{HPD}{logical value. If \code{TRUE} then only those subjects are
    classified for which the lower limit of the \code{pHPD}*100\% HPD
    credible interval of the component probability exceeds the value of \code{ptrash}.}
  \item{pHPD}{credible level of the HPD credible interval, see argument \code{HPD}.}
  \item{pthresh}{an optional threshold for the estimated component
    probability (when \code{HPD} is \code{FALSE}) or for the lower limit
    of the HPD credible interval (when \code{HPD} is \code{TRUE}) to
    classify a subject. No effect when \code{pthresh} is negative.}
  \item{unclass.na}{logical value taken into account when \code{pthresh}
   is positive. If \code{unclass.na} is \code{TRUE}, unclassified subjects get classification
   \code{NA}. If \code{unclass.na} is \code{FALSE}, unclassified
   subjects create a separate (last) group.}
  \item{\dots}{optional additional arguments.}  
}  
\value{
  A \code{data.frame} with three (when \code{HPD} is \code{FALSE}) or five
  (when \code{HPD} is \code{TRUE}) columns.
}  
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
## TO BE ADDED.
}
\keyword{multivariate}
\keyword{cluster}
