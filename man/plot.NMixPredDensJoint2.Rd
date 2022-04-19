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
    add.contour=TRUE, col.add.contour="brown",
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
  \item{add.contour}{logical. If \code{TRUE} and \code{contour} is
    \code{FALSE} (i.e., image plot is drawn) then contours are added to
    the image plots.}
  \item{col.add.contour}{color of contours which are added to the image plot.}
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
  Komárek, A. (2009).
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis}, \bold{53}(12), 3932--3947.
}
\seealso{
  \code{\link{NMixPredDensJoint2}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
## See additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## - files Galaxy.R, Faithful.R, Tandmob.R and
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Galaxy.pdf
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Faithful.pdf
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Tandmob.pdf
}
\keyword{dplot}
