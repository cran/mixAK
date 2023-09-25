\name{NMixPredCondDensMarg}
\alias{NMixPredCondDensMarg}
\alias{NMixPredCondDensMarg.default}
\alias{NMixPredCondDensMarg.NMixMCMC}
\alias{NMixPredCondDensMarg.GLMM_MCMC}
\alias{C_NMix_PredCondDensCDFMarg}
\title{
  Univariate conditional predictive density
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  (posterior predictive) estimates of univariate conditional densities.
}
\usage{
NMixPredCondDensMarg(x, \dots)

\method{NMixPredCondDensMarg}{default}(x, icond, prob, scale, K, w, mu, Li, Krandom=FALSE, ...)

\method{NMixPredCondDensMarg}{NMixMCMC}(x, icond, prob, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPredCondDensMarg}{GLMM_MCMC}(x, icond, prob, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredCondDensMarg.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPredCondDensMarg.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPredCondDensMarg.default} function.
  }
  \item{icond}{index of the margin by which we want to condition
  }
  \item{prob}{a numeric vector. If given then also the posterior
    pointwise quantiles of the conditional densities are computed for
    probabilities given by \code{prob}. These can be used to draw
    pointwise credible intervals.
  }  
  \item{scale}{a two component list giving the \code{shift} and the
    \code{scale}. If not given, shift is equal to zero and scale is
    equal to one.
  }
  \item{K}{either a number (when \code{Krandom}\eqn{=}\code{FALSE}) or a
    numeric vector with the chain for the number of mixture components.
  }
  \item{w}{a numeric vector with the chain for the mixture weights.
  }
  \item{mu}{a numeric vector with the chain for the mixture means.
  }
  \item{Li}{a numeric vector with the chain for the mixture inverse
    variances (lower triangles only).
  }
  \item{Krandom}{a logical value which indicates whether the number of
    mixture components changes from one iteration to another.
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
  An object of class \code{NMixPredCondDensMarg} which has the following components:  
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
  \item{prob}{a value of the argument \code{prob}.
  }
  \item{qXX\%}{if \code{prob} is given then there is one
    additional component named qXX\%, e.g., q50\% for
    each value of \code{prob} which has the same structure as the
    component \code{dens} and keeps computed posterior pointwise
    quantiles.
  }  
  There is also a \code{plot} method implemented for the resulting object.
}
\seealso{
  \code{\link{plot.NMixPredCondDensMarg}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
