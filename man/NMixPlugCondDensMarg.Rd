\name{NMixPlugCondDensMarg}
\alias{NMixPlugCondDensMarg}
\alias{NMixPlugCondDensMarg.default}
\alias{NMixPlugCondDensMarg.NMixMCMC}
\alias{NMixPlugCondDensMarg.GLMM_MCMC}
\title{
  Univariate conditional densities: plug-in estimate
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  estimates of univariate conditional densities obtained by using posterior
  summary statistics (e.g., posterior means) of mixture weights, means
  and variances (plug-in estimate).
}
\usage{
NMixPlugCondDensMarg(x, \dots)

\method{NMixPlugCondDensMarg}{default}(x, icond, scale, w, mu, Sigma, \dots)

\method{NMixPlugCondDensMarg}{NMixMCMC}(x, icond, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPlugCondDensMarg}{GLMM_MCMC}(x, icond, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPlugCondDensMarg.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPlugCondDensMarg.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPlugCondDensMarg.default} function.
  }
  \item{icond}{index of the margin by which we want to condition
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
    the density should be evaluated. The value of \code{grid[[icond]]}
    determines the values by which we condition.

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
  An object of class \code{NMixPlugCondDensMarg} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{icond}{index of the margin by which we condition.}
  \item{dens}{a list with the computed conditional densities for each
    value of \code{x[[icond]]}. Each \code{dens[[j]]} is again a list
    with conditional densities for each margin given margin
    \code{icond} equal to \code{x[[icond]][j]}.
    The value of \code{dens[[j]][[imargin]]} gives a value
    of a marginal density of the \code{imargin}-th margin at \code{x[[icond]][j]}.
  }    
  There is also a \code{plot} method implemented for the resulting object.
}
\seealso{
  \code{\link{plot.NMixPlugCondDensMarg}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}, \code{\link{NMixPredCondDensMarg}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
