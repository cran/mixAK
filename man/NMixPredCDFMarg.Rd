\name{NMixPredCDFMarg}
\alias{NMixPredCDFMarg}
\alias{NMixPredCDFMarg.default}
\alias{NMixPredCDFMarg.NMixMCMC}
\alias{NMixPredCDFMarg.GLMM_MCMC}
\title{
  Marginal (univariate) predictive cumulative distribution function
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  estimated posterior predictive cumulative distribution function for each margin.
}
\usage{
NMixPredCDFMarg(x, \dots)

\method{NMixPredCDFMarg}{default}(x, scale, K, w, mu, Li, Krandom=TRUE, \dots)

\method{NMixPredCDFMarg}{NMixMCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPredCDFMarg}{GLMM_MCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredCDFMarg.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPredCDFMarg.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPredCDFMarg.default} function.
  }
  \item{scale}{a two component list giving the \code{shift} and the \code{scale}.}
  \item{K}{either a number (when \code{Krandom}\eqn{=}\code{FALSE}) or a
    numeric vector with the chain for the number of mixture components.}
  \item{w}{a numeric vector with the chain for the mixture weights.}
  \item{mu}{a numeric vector with the chain for the mixture means.}
  \item{Li}{a numeric vector with the chain for the mixture inverse
    variances (lower triangles only).}
  \item{Krandom}{a logical value which indicates whether the number of
    mixture components changes from one iteration to another.}
  \item{grid}{a numeric vector or a list with the grid values in which
    the predictive CDF should be evaluated.

    If \code{x$dim} is 1 then \code{grid} may be a numeric vector. If
    \code{x$dim} is higher than then \code{grid} must be a list with
    numeric vectors as components giving the grids for each margin.

    If \code{grid} is not specified, it is created automatically using
    the information from the posterior summary statistics stored in \code{x}.
  }
  \item{lgrid}{a length of the grid used to create the \code{grid} if
    that is not specified.}
  \item{scaled}{if \code{TRUE}, the CDF of shifted and scaled data is
    summarized. The shift and scale vector are taken from the
    \code{scale} component of the object \code{x}.}
  \item{\dots}{optional additional arguments.}
}
\value{
  An object of class \code{NMixPredCDFMarg} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{freqK}{frequency table for the values of \eqn{K} (numbers of
    mixture components) in the MCMC chain.}
  \item{propK}{proportions derived from \code{freqK}.}
  \item{MCMC.length}{the length of the MCMC used to compute the
    predictive cdf's.}
  \item{cdf}{a list with the computed predictive CDF's for each
    margin. The components of the list are named \code{1}, \ldots, i.e.,
    \code{cdf[[1]]}\eqn{=}\code{cdf[["1"]]} is the predictive
    cdf for margin 1 etc.}
  \item{cdfK}{a list with the computed predictive CDF's for each
    margin, conditioned further by \eqn{K}. The components of the list
    are named \code{1}, \ldots. That is,
    \code{cdf[[1]][[1]]} \eqn{=} \code{cdf[["1"]][[1]]} is the predictive
    CDF for margin 1 conditioned by \eqn{K=1}{K=1},
    \code{cdf[[1]][[2]]} \eqn{=} \code{cdf[["1"]][[2]]} is the predictive
    CDF for margin 1 conditioned by \eqn{K=2}{K=2} etc.

    Note that \code{cdfK} provides some additional information only
    when \code{Krandom} \eqn{=} \code{TRUE} or when \code{x} results from
    the \code{\link{NMixMCMC}} call to the reversible jump MCMC.
  }
  There is also a \code{plot} method implemented for the resulting object.
}
\references{
  Komárek, A. (2009).  
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis}, \bold{53}(12), 3932--3947.
}
\seealso{
  \code{\link{plot.NMixPredCDFMarg}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
