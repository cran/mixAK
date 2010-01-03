\name{NMixPredDensJoint2}
\alias{NMixPredDensJoint2}
\alias{NMixPredDensJoint2.default}
\alias{NMixPredDensJoint2.NMixMCMC}
\alias{NMixPredDensJoint2.GLMM_MCMC}
\title{
  Pairwise bivariate predictive density
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  estimated posterior predictive densities for each pair of margins.
}
\usage{
NMixPredDensJoint2(x, \dots)

\method{NMixPredDensJoint2}{default}(x, scale, K, w, mu, Li, Krandom=TRUE, \dots)

\method{NMixPredDensJoint2}{NMixMCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPredDensJoint2}{NMixMCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredDensJoint2.NMixMCMC} function.

    an object of class \code{NMixMCMC} for
    \code{NMixPlugDensJoint2.NMixMCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPredDensJoint2.default} function.
  }
  \item{scale}{a two component list giving the \code{shift} and the
    \code{scale}. If not given, shift is equal to zero and scale is
    equal to one.}
  \item{K}{either a number (when \code{Krandom}\eqn{=}\code{FALSE}) or a
    numeric vector with the chain for the number of mixture components.}
  \item{w}{a numeric vector with the chain for the mixture weights.}
  \item{mu}{a numeric vector with the chain for the mixture means.}
  \item{Li}{a numeric vector with the chain for the mixture inverse
    variances (lower triangles only).}
  \item{Krandom}{a logical value which indicates whether the number of
    mixture components changes from one iteration to another.}
  \item{grid}{a list with the grid values for each margin in which
    the predictive density should be evaluated.

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
  An object of class \code{NMixPredDensJoint2} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{freqK}{frequency table for the values of \eqn{K} (numbers of
    mixture components) in the MCMC chain.}
  \item{propK}{proportions derived from \code{freqK}.}
  \item{MCMC.length}{the length of the MCMC used to compute the
    predictive densities.}
  \item{dens}{a list with the computed predictive densities for each
    pair of margins. The components of the list are named \code{1-2},
    \code{1-3}, \ldots, i.e.,
    \code{dens[[1]]}\eqn{=}\code{dens[["1-2"]]} is the pairwise predictive
    density for margins 1 and 2, etc. Each component of the \code{list}
    is a matrix in such a form that it can be directly passed together
    with the proper components of \code{x} to the plotting functions
    like \code{\link[graphics]{contour}} or \code{\link[graphics]{image}}.
  }
  \item{densK}{a list with the computed predictive densities for each
    margin, conditioned further by \eqn{K}. The components of the list
    are named \code{1-2}, \code{1-3}, \ldots. That is,
    \code{dens[[1]][[1]]} \eqn{=} \code{dens[["1-2"]][[1]]} is the
    pairwise predictive  density for margins 1 and 2 conditioned by \eqn{K=1}{K=1},
    \code{dens[[1]][[2]]} \eqn{=} \code{dens[["1-2"]][[2]]} is the
    pairwise predictive density for margins 1 and 2 conditioned by \eqn{K=2}{K=2} etc.

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
  \emph{Computational Statistics and Data Analysis}, \bold{53}, 3932--3947.
}
\seealso{
  \code{\link{plot.NMixPredDensJoint2}}, \code{\link{NMixMCMC}},
  \code{\link{GLMM_MCMC}}, \code{\link{NMixPredDensMarg}}.
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


