\name{NMixPredCondCDFMarg}
\alias{NMixPredCondCDFMarg}
\alias{NMixPredCondCDFMarg.default}
\alias{NMixPredCondCDFMarg.NMixMCMC}
\alias{NMixPredCondCDFMarg.GLMM_MCMC}
\title{
  Univariate conditional predictive cumulative distribution function
}
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  (posterior predictive) estimates of univariate conditional cumulative distribution functions.
}
\usage{
NMixPredCondCDFMarg(x, \dots)

\method{NMixPredCondCDFMarg}{default}(x, icond, prob, scale, K, w, mu, Li, Krandom=FALSE, ...)

\method{NMixPredCondCDFMarg}{NMixMCMC}(x, icond, prob, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPredCondCDFMarg}{GLMM_MCMC}(x, icond, prob, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredCondCDFMarg.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPredCondCDFMarg.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPredCondCDFMarg.default} function.
  }
  \item{icond}{index of the margin by which we want to condition
  }
  \item{prob}{a numeric vector. If given then also the posterior
    pointwise quantiles of the conditional cdf's are computed for
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
    the cdf should be evaluated. The value of \code{grid[[icond]]}
    determines the values by which we condition.

    If \code{grid} is not specified, it is created automatically using
    the information from the posterior summary statistics stored in \code{x}.
  }
  \item{lgrid}{a length of the grid used to create the \code{grid} if
    that is not specified.
  }
  \item{scaled}{if \code{TRUE}, the cdf of shifted and scaled data is
    summarized. The shift and scale vector are taken from the
    \code{scale} component of the object \code{x}.
  }
  \item{\dots}{optional additional arguments.}    
}
\value{
  An object of class \code{NMixPredCondCDFMarg} which has the following components:  
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{icond}{index of the margin by which we condition.}
  \item{cdf}{a list with the computed conditional cdf's for each
    value of \code{x[[icond]]}. Each \code{cdf[[j]]} is again a list
    with conditional cdf's for each margin given margin
    \code{icond} equal to \code{x[[icond]][j]}.
    The value of \code{cdf[[j]][[imargin]]} gives a value
    of a marginal cdf of the \code{imargin}-th margin at \code{x[[icond]][j]}.
  }
  \item{prob}{a value of the argument \code{prob}.
  }
  \item{qXX\%}{if \code{prob} is given then there is one
    additional component named qXX\%, e.g., q50\% for
    each value of \code{prob} which has the same structure as the
    component \code{cdf} and keeps computed posterior pointwise
    quantiles.
  }
  There is also a \code{plot} method implemented for the resulting object.
}
\seealso{
  \code{\link{plot.NMixPredCondCDFMarg}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
