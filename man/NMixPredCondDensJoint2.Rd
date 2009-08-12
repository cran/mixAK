\name{NMixPredCondDensJoint2}
\alias{NMixPredCondDensJoint2}
\alias{NMixPredCondDensJoint2.default}
\alias{NMixPredCondDensJoint2.NMixMCMC}
\alias{NMixPredCondDensJoint2.GLMM_MCMC}
\title{
  Pairwise bivariate conditional predictive densities
}  
\description{
  This function serves as an inference tool for the MCMC output
  obtained using the function \code{\link{NMixMCMC}}. It computes
  estimates of pairwise bivariate conditional densities (given one margin)
  obtained by using posterior summary statistics (e.g., posterior means)
  of mixture weights, means and variances (plug-in estimate).
}
\usage{
NMixPredCondDensJoint2(x, \dots)

\method{NMixPredCondDensJoint2}{default}(x, icond, scale, K, w, mu, Li, Krandom=FALSE, \dots)

\method{NMixPredCondDensJoint2}{NMixMCMC}(x, icond, grid, lgrid=50, scaled=FALSE, \dots)

\method{NMixPredCondDensJoint2}{GLMM_MCMC}(x, icond, grid, lgrid=50, scaled=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} for
    \code{NMixPredCondDensJoint2.NMixMCMC} function.

    An object of class \code{GLMM_MCMC} for
    \code{NMixPredCondDensJoint2.GLMM_MCMC} function.
    
    A list with the grid values (see below) for
    \code{NMixPredCondDensJoint2.default} function.
  }
  \item{icond}{index of the margin by which we want to condition
  } 
  \item{scale}{a two component list giving the \code{shift} and the
    \code{scale}. If not given, shift is equal to zero and scale is
    equal to one.
  }
  \item{K}{either a number (when \code{Krandom}\eqn{=}\code{FALSE}) or a
    numeric vector with the chain for the number of mixture components.}
  \item{w}{a numeric vector with the chain for the mixture weights.}
  \item{mu}{a numeric vector with the chain for the mixture means.}
  \item{Li}{a numeric vector with the chain for the mixture inverse
    variances (lower triangles only).}
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
  An object of class \code{NMixPredCondDensJoint2} which has the
  following components:
  \item{x}{a list with the grid values for each margin. The components
    of the list are named \code{x1}, \ldots or take names from
    \code{grid} argument.}
  \item{icond}{index of the margin by which we condition.}
  \item{dens}{a list with the computed conditional densities for each
    value of \code{x[[icond]]}. Each \code{dens[[j]]} is again a list
    with conditional densities for each pair of margins given margin
    \code{icond} equal to \code{x[[icond]][j]}.
    The value of \code{dens[[j]][[\dQuote{i-k}]]} gives values
    of conditional density of the (i,k)-th margins given margin
    \code{icond} equal to \code{x[[icond]][j]}.
  }    
  There is also a \code{plot} method implemented for the resulting object.  
}
\seealso{
  \code{\link{plot.NMixPredCondDensJoint2}}, \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{dplot}
\keyword{smooth}
