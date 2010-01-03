\name{GLMM_longitDA}
\alias{GLMM_longitDA}
\title{
  Discriminant analysis for longitudinal profiles based on fitted GLMM's
}  
\description{
  THIS FUNCTION IS BEING DEVELOPED AND ORDINARY USERS
  ARE NOT RECOMMENDED TO USE IT.
}
\details{
  I HOPE TO WRITE A PAPER DESCRIBING ALL THE DETAILS.
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
  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.}, Hansen, B. E., Kuiper,
  E. M. M., van Buuren, H. R., and Lesaffre, E. (2010).
  Discriminant analysis using a multivariate linear mixed model with a
  normal mixture in the random effects distribution.
  \emph{Statistics in Medicine}. To appear. %\bold{XX}, XXX-XXX.
}  
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### WILL BE ADDED.
}
\keyword{models}
\keyword{multivariate}
\keyword{cluster}
