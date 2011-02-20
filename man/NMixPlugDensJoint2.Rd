\name{NMixPlugDensJoint2}
\alias{NMixPlugDensJoint2}
\alias{NMixPlugDensJoint2.default}
\alias{NMixPlugDensJoint2.NMixMCMC}
\alias{NMixPlugDensJoint2.GLMM_MCMC}
\title{
  Pairwise bivariate densities: plug-in estimate
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  marginal (pairwise bivariate) plug-in densities obtained by using posterior
  summary statistics (e.g., posterior means) of mixture weights, means
  and variances.
}
\usage{
NMixPlugDensJoint2(x, \dots)

\method{NMixPlugDensJoint2}{default}(x, scale, w, mu, Sigma, \dots)

\method{NMixPlugDensJoint2}{NMixMCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPlugDensJoint2}{GLMM_MCMC}(x, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPlugDensJoint2.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPlugDensJoint2.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPlugDensJoint2.default} function.
  }
  \item{scale}{a two component list giving the \code{shift} and the
    \code{scale}. If not given, shift is equal to zero and scale is
    equal to one.
  }
  \item{w}{a numeric vector with posterior summary statistics for the
    mixture weights. The length of this vector determines the number of
    mixture components.
  }
  \item{mu}{a matrix with posterior summary statistics for
    mixture means in rows. That is, \code{mu} has
    \eqn{K} rows and \eqn{p} columns, where \eqn{K} denotes the number
    of mixture components and \eqn{p} is dimension of the mixture
    distribution.
  }
  \item{Sigma}{a list with posterior summary statistics for for mixture covariance matrices.
  }
  \item{grid}{a list with the grid values for each margin in which
    the density should be evaluated.

    If \code{grid} is not specified, it is created automatically using
    the information from the posterior summary statistics stored in \code{x}.
  }
  \item{lgrid}{a length of the grid used to create the \code{grid} if
    that is not specified.
  }
  \item{scaled}{if \code{TRUE}, the density of shifted and scaled data is
    summarized. The shift and scale vector are taken from the
    \code{scale} component of the object \code{x}.
  }
  \item{\dots}{optional additional arguments.}    
}
\value{
  An object of class \code{NMixPlugDensJoint2} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{dens}{a list with the computed densities for each
    pair of margins. The components of the list are named \code{1-2},
    \code{1-3}, \ldots, i.e.,
    \code{dens[[1]]}\eqn{=}\code{dens[["1-2"]]} is the pairwise predictive
    density for margins 1 and 2, etc. Each component of the \code{list}
    is a matrix in such a form that it can be directly passed together
    with the proper components of \code{x} to the plotting functions
    like \code{\link[graphics]{contour}} or \code{\link[graphics]{image}}.
  }
  There is also a \code{plot} method implemented for the resulting object.
}
\seealso{
  \code{\link{plot.NMixPlugDensJoint2}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}, \code{\link{NMixPredDensJoint2}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
