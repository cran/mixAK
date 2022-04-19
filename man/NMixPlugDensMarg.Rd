\name{NMixPlugDensMarg}
\alias{NMixPlugDensMarg}
\alias{NMixPlugDensMarg.default}
\alias{NMixPlugDensMarg.NMixMCMC}
\alias{NMixPlugDensMarg.GLMM_MCMC}
\title{
  Marginal (univariate) densities: plug-in estimate
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  marginal (univariate) plug-in densities obtained by using posterior
  summary statistics (e.g., posterior means) of mixture weights, means
  and variances.
}
\usage{
NMixPlugDensMarg(x, \dots)

\method{NMixPlugDensMarg}{default}(x, scale, w, mu, Sigma, \dots)

\method{NMixPlugDensMarg}{NMixMCMC}(x, grid, lgrid=500, scaled=FALSE, \dots)

\method{NMixPlugDensMarg}{GLMM_MCMC}(x, grid, lgrid=500, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{a list with the grid values (see below) for
    \code{NMixPlugDensMarg.default} function.

    An object of class \code{NMixMCMC} for
    \code{NMixPlugDensMarg.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPlugDensMarg.GLMM_MCMC} function.    
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
  \item{Sigma}{a list with posterior summary statistics for mixture covariance matrices.
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
  An object of class \code{NMixPlugDensMarg} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{dens}{a list with the computed densities for each
    margin. The components of the list are named \code{1}, \ldots, i.e.,
    \code{dens[[1]]}\eqn{=}\code{dens[["1"]]} is the predictive
    density for margin 1 etc.}
  There is also a \code{plot} method implemented for the resulting object.
}
\seealso{
  \code{\link{plot.NMixPlugDensMarg}}, \code{\link{NMixMCMC}},
  \code{\link{GLMM_MCMC}}, \code{\link{NMixPredDensMarg}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
