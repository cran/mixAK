\name{NMixPredDensMarg}
\alias{NMixPredDensMarg}
\alias{NMixPredDensMarg.default}
\alias{NMixPredDensMarg.NMixMCMC}
\title{
  Marginal (univariate) predictive density
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  estimated posterior predictive densities for each margin.
}
\usage{
NMixPredDensMarg(x, \dots)

\method{NMixPredDensMarg}{default}(x, scale, K, w, mu, Li, Krandom=TRUE, \dots)

\method{NMixPredDensMarg}{NMixMCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredDensMarg.NMixMCMC} function.

    A list with the grid values (see below) for
    \code{NMixPredDensMarg.default} function.
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
    the predictive density should be evaluated.

    If \code{x$dim} is 1 then \code{grid} may be a numeric vector. If
    \code{x$dim} is higher than then \code{grid} must be a list with
    numeric vectors as components giving the grids for each margin.

    If \code{grid} is not specified, it is created automatically using
    the information from the posterior summary statistics stored in \code{x}.
  }
  \item{lgrid}{a length of the grid used to create the \code{grid} if
    that is not specified.}
  \item{scaled}{if \code{TRUE}, the density of shifted and scaled data is
    summarized. The shift and scale vector are taken from the
    \code{scale} component of the object \code{x}.}
  \item{\dots}{optional additional arguments.}
}
\value{
  An object of class \code{NMixPredDensMarg} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots.}
  \item{freqK}{frequency table for the values of \eqn{K} (numbers of
    mixture components) in the MCMC chain.}
  \item{propK}{proportions derived from \code{freqK}.}
  \item{MCMC.length}{the length of the MCMC used to compute the
    predictive densities.}
  \item{dens}{a list with the computed predictive densities for each
    margin. The components of the list are named \code{1}, \ldots, i.e.,
    \code{dens[[1]]}\eqn{=}\code{dens[["1"]]} is the predictive
    density for margin 1 etc.}
  \item{densK}{a list with the computed predictive densities for each
    margin, conditioned further by \eqn{K}. The components of the list
    are named \code{1}, \ldots. That is,
    \code{dens[[1]][[1]]} \eqn{=} \code{dens[["1"]][[1]]} is the predictive
    density for margin 1 conditioned by \eqn{K=1}{K=1},
    \code{dens[[1]][[2]]} \eqn{=} \code{dens[["1"]][[2]]} is the predictive
    density for margin 1 conditioned by \eqn{K=2}{K=2} etc.

    Note that \code{densK} provides some additional information only
    when \code{Krandom} \eqn{=} \code{TRUE} or when \code{x} results from
    the \code{\link{NMixMCMC}} call to the reversible jump MCMC.
  }
  There is also a \code{plot} method implemented for the resulting object.
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.}
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis.} To appear.  
}
\seealso{
  \code{\link{plot.NMixPredDensMarg}}, \code{\link{NMixMCMC}}, \code{\link{NMixPredDensJoint2}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
## See additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## - files Galaxy.pdf, Faithful.pdf, Tandmob.pdf
}
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
