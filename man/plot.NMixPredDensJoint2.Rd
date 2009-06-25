\name{plot.NMixPredDensJoint2}
\alias{plot.NMixPredDensJoint2}
\title{
  Plot computed marginal pairwise bivariate predictive densities
}
\description{
  This is a basic plotting tool to visualize computed marginal pairwise
  bivariate predictive densities using the
  \code{\link[graphics]{image}} plot or \code{\link[graphics]{contour}}
  plot.
  See also \code{\link{NMixPredDensJoint2}}.
}
\usage{
\method{plot}{NMixPredDensJoint2}(x, K=0, contour=FALSE,
    auto.layout=TRUE,
    col, lwd=1, main, xylab, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPredDensJoint2}.}
  \item{K}{if equal to \code{0} then the overall predictive densities
    are plotted taken from the \code{dens} part of the object \code{x}.

    If higher than \code{0} then the predictive density conditioned by
    the value of \code{K} is plotted (taken from the \code{densK} part
    of the object \code{x}).
  }
  \item{contour}{logical. If \code{TRUE} then contours are drawn,
    otherwise image plot is created.
  }
  \item{auto.layout}{if \code{TRUE} then the function determines itself
    how to divide the plotting region to draw densities for all margins.
  }
  \item{col}{color used to draw the contours or images.}
  \item{lwd}{line width.}
  \item{main}{main title of the plot.}
  \item{xylab}{optional character vector of the length equal to
    the number of margins with labels used for x and y axes on the plots.}
  \item{\dots}{additional arguments passed to the \code{plot} function.} 
}
\value{
  \code{invisible(x)}
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.}
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis.} To appear.  
}
\seealso{
  \code{\link{NMixPredDensJoint2}}, \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
## See additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## - files Galaxy.pdf, Faithful.pdf, Tandmob.pdf
}
\keyword{dplot}
