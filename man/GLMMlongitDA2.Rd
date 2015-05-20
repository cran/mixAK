\name{GLMM_longitDA2}
\alias{GLMM_longitDA2}
\title{
  Discriminant analysis for longitudinal profiles based on fitted GLMM's
}  
\description{
  WILL BE ADDED.
}
\details{
  This function complements a paper being currently in preparation.

  \code{GLMM_longitDA2} differs in many aspects from \code{\link{GLMM_longitDA2}}!
}  
\usage{
GLMM_longitDA2(mod, w.prior, y, id, x, z, xz.common = TRUE,
               keep.comp.prob = FALSE, level = 0.95,
               info, silent = FALSE)
}
\arguments{
  \item{mod}{a list containing models fitted with the
    \code{\link{GLMM_MCMC}} function. Each component of the list is the
    GLMM fitted in the training dataset of each cluster.
  }
  \item{w.prior}{a vector with prior cluster weights. The length of this
    argument must be the same as the length of argument \code{mod}.
    Can also be given relatively, e.g., as \code{c(1, 1)} which means
    that both prior weights are equal to 1/2.
  }
  \item{y}{vector, matrix or data frame (see argument \code{y} of
    \code{\link{GLMM_MCMC}} function) with responses of objects that are
    to be clustered.
  }
  \item{id}{vector which determines clustered observations (see also
    argument \code{y} of \code{\link{GLMM_MCMC}} function).
  }
  \item{x}{see \code{xz.common} below.}
  \item{z}{see \code{xz.common} below.}
  \item{xz.common}{a logical value.

    If \code{TRUE} then it is assumed
    that the X and Z matrices are the same for GLMM in each cluster. In
    that case, arguments \code{x} and \code{z} have the same structure
    as arguments \code{x} and \code{z} of \code{\link{GLMM_MCMC}}
    function.

    If \code{FALSE} then X and Z matrices for the GLMM may differ across
    clusters. In that case, arguments \code{x} and \code{z} are both
    lists of length equal to the number of clusters and each component
    of lists \code{x} and \code{z} has the same structure as arguments
    \code{x} and \code{z} of \code{\link{GLMM_MCMC}} function.
  }
  \item{keep.comp.prob}{a logical value indicating whether the
    allocation probabilities should be kept for all MCMC
    iterations. This may ask for quite some memory but it is necessary
    if credible intervals etc. should be calculated for the component
    probabilities.}
  \item{level}{level of HPD credible intervals that are calculated for
    the component probabilities if \code{keep.comp.prob} is \code{TRUE}.}
  \item{info}{interval in which the function prints the progress of
    computation (unless \code{silent} is \code{TRUE}).}
  \item{silent}{logical value indicating whether to switch-off printing
    the information during calculations.}
}
\value{
  A list with the following components:
  \item{ADD}{ADD DESCRIPTION}
}
\seealso{
  \code{\link{GLMM_MCMC}}, \code{\link{GLMM_longitDA}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### WILL BE ADDED.
}
\keyword{models}
\keyword{multivariate}
\keyword{cluster}
