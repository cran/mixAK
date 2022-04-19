\name{plot.NMixPlugDensJoint2}
\alias{plot.NMixPlugDensJoint2}
\title{
  Plot computed marginal pairwise bivariate densities (plug-in estimate)
}
\description{
  This is a basic plotting tool to visualize computed marginal pairwise
  bivariate densities (plug-in version) using the
  \code{\link[graphics]{contour}} plot. See also \code{\link{NMixPlugDensJoint2}}.
}
\usage{
\method{plot}{NMixPlugDensJoint2}(x, contour=FALSE,
    add.contour=TRUE, col.add.contour="brown",
    auto.layout=TRUE,
    col, lwd=1, main, xylab, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPlugDensJoint2}.}
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
    the number of margins with labels used for x and y axes on the
    plots.
  }
  \item{\dots}{additional arguments passed to the \code{plot} function.} 
}
\value{
  \code{invisible(x)}
}
\seealso{
  \code{\link{NMixPlugDensJoint2}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{dplot}
