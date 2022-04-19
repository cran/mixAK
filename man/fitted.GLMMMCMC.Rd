\name{fitted.GLMM_MCMC}
\alias{fitted.GLMM_MCMC}
\title{Fitted profiles in the GLMM model}
\description{It calculates fitted profiles in the (multivariate) GLMM
  with a normal mixture in the random effects distribution based on
  selected posterior summary statistic of the model parameters.}
\usage{
\method{fitted}{GLMM_MCMC}(object, x, z,
  statistic=c("median", "mean", "Q1", "Q3", "2.5\%", "97.5\%"),
  overall=FALSE, glmer=TRUE, nAGQ=100, \dots)
}
\arguments{
  \item{object}{object of class \code{\link{GLMM_MCMC}}.}
  \item{x}{matrix or list of matrices (in the case of multiple
    responses) for \dQuote{fixed effects} part of the model used in the
    calculation of fitted values.}
  \item{z}{matrix or list of matrices (in the case of multiple
    responses) for \dQuote{random effects} part of the model used in the
    calculation of fitted values.}
  \item{statistic}{character which specifies the posterior summary
    statistic to be used to calculate fitted profiles. Default is the
    posterior median. It applies only to the overall fit.}
  \item{overall}{logical. If \code{TRUE}, fitted profiles based on posterior
    mean/median/Q1/Q3/2.5\%/97.5\% of the model parameters are computed.
    
    If \code{FALSE}, fitted profiles based on posterior means
    given mixture component are calculated. Note that this depends on
    used re-labelling of the mixture components and hence might be
    misleading if re-labelling is not succesfull!}
  \item{glmer}{a logical value. If \code{TRUE}, the real marginal means
    are calculated using Gaussian quadrature.}
  \item{nAGQ}{number of quadrature points used when \code{glmer} is
    \code{TRUE}.}
  \item{\dots}{possibly extra arguments. Nothing useful at this moment.}
}
\value{
  A list (one component for each of multivariate responses from the
  model) with fitted values calculated using the \code{x} and \code{z}
  matrices. If \code{overall} is \code{FALSE}, these are then matrices
  with one column for each mixture component.
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
### WILL BE ADDED.
}
\keyword{smooth}
