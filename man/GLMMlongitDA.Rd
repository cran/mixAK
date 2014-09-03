\name{GLMM_longitDA}
\alias{GLMM_longitDA}
\title{
  Discriminant analysis for longitudinal profiles based on fitted GLMM's
}  
\description{
  The idea is that we fit (possibly different) GLMM's for data in training 
  groups using the function \code{\link{GLMM_MCMC}} and then use the fitted 
  models for discrimination of new observations. For more details we 
  refer to Komárek et al. (2010).
}
\details{
  This function complements a paper Komárek et al. (2010).
}  
\usage{
GLMM_longitDA(mod, w.prior, y, id, time, x, z, xz.common=TRUE, info)
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
  \item{time}{vector which gives indeces of observations within
    clusters. It appears (together with \code{id}) in the output as
    identifier of observations}  
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
  \item{info}{interval in which the function prints the progress of computation}
}
\value{
  A list with the following components:
  \item{ident}{ADD DESCRIPTION}
  \item{marg}{ADD DESCRIPTION}
  \item{cond}{ADD DESCRIPTION}
  \item{ranef}{ADD DESCRIPTION}
}
\references{
  Komárek, A., Hansen, B. E., Kuiper,
  E. M. M., van Buuren, H. R., and Lesaffre, E. (2010).
  Discriminant analysis using a multivariate linear mixed model with a
  normal mixture in the random effects distribution.
  \emph{Statistics in Medicine}, \bold{29}(30), 3267--3283.
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
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
