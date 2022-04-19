\name{NMixPseudoGOF}
\alias{NMixPseudoGOF}
\alias{NMixPseudoGOF.default}
\alias{NMixPseudoGOF.NMixMCMC}
\title{
  Pseudo goodness-of-fit test for a normal mixture model
}  
\description{
  It takes a (fitted) normal mixture, creates hyperrectangles according to a specified grid,
  computes probability masses in each hyperrectangle derived from the (fitted) normal mixture.
  From computed probability masses expected frequencies (using the sample size of supplied data) are
  computed and compared to frequencies observed in supplied data. From expected and observed
  frequencies, a Pearson chi-squared like statistic is computed and returned together with
  residuals derived from that statistic.

  Also pseudo degrees of freedom are returned which are equal
  to a number of hyperrectangles minus number of free parameters of the normal mixture.
  For a \eqn{K}-component mixture of dimension \eqn{p}, the number of free parameters
  is computed as
  \deqn{q = K-1 + K\cdot p + K\cdot p(p+1)/2}{q = K-1 + K*p + K*p*(p+1)/2}
  Note that computation of \eqn{q} does not take into account the positive (semi-)definiteness
  restriction on covariance matrices.

  WARNING: There is no statistical theory developed that would guarantee that
  computed chi-squared like statistics follows a chi-squared distribution
  with computed pseudo degrees of freedom under the null hypothesis that
  the distribution that generated the data is a normal mixture. This function serves purely
  for descriptive purposes!  
}  
\usage{
NMixPseudoGOF(x, \dots)

\method{NMixPseudoGOF}{default}(x, scale, w, mu, Sigma, breaks, nbreaks=10, digits=3, \dots)

\method{NMixPseudoGOF}{NMixMCMC}(x, y, breaks, nbreaks=10, digits=3, \dots)
}
\arguments{
  \item{x}{data object (see argument \code{y} below) for
    \code{NMixPseudoGOF.default} function.

    An object of class \code{NMixMCMC} for
    \code{NMixPseudoGOF.NMixMCMC} function.
  }
  \item{y}{a numeric vector, matrix or data frame with the data. It is a
    numeric vector if \eqn{p} is one. It is a matrix or data frame with
    \eqn{p} columns if \eqn{p > 1}{p > 1}.
  }
  \item{scale}{a two component list giving the \code{shift} and the
    \code{scale}. If not given, shift is equal to zero and scale is
    equal to one.
  }
  \item{w}{a numeric vector with mixture weights.
    The length of this vector determines the number of
    mixture components.
  }
  \item{mu}{a matrix with mixture means in rows. That is, \code{mu} has
    \eqn{K} rows and \eqn{p} columns, where \eqn{K} denotes the number
    of mixture components and \eqn{p} is dimension of the mixture
    distribution.
  }
  \item{Sigma}{a list with mixture covariance matrices.
  }
  \item{breaks}{a numeric vector or a list with the breaks defining the
    hyperrectangles. It is a numeric vector if \eqn{p} is equal to
    one. It is a list of length \eqn{p} of numeric vectors. Each
    component of the list determines the breaks for each margin.     
  }
  \item{nbreaks}{a number or a numeric vector with the number of breaks
    for each margin. It is only used if the argument \code{breaks} is
    not given to determine sensible break values.
  }
  \item{digits}{a number or a numeric vector with the number of digits
    to which the breaks should be rounded in the case they are created
    by the function. If it is a vector then different rounding may be
    used for each margin.
  }
  \item{\dots}{optional additional arguments.}  
}
\value{
ADD DESCRIPTION
}
\seealso{
  \code{\link{NMixMCMC}}.    
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
