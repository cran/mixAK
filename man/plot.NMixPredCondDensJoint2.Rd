\name{plot.NMixPredCondDensJoint2}
\alias{plot.NMixPredCondDensJoint2}
\title{
  Plot computed predictive pairwise bivariate conditional densities
}
\description{
  This is a basic plotting tool to visualize computed 
  predictive pairwise bivariate conditional densities using the
  \code{\link[graphics]{image}} or \code{\link[graphics]{contour}} plot.
  See also \code{\link{NMixPredCondDensJoint2}}.
}
\usage{
\method{plot}{NMixPredCondDensJoint2}(x, ixcond, imargin,
    contour=FALSE, auto.layout=TRUE,
    col, lwd=1, main, xylab, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPredCondDensJoint2}.}
  \item{ixcond}{if given then conditional densities of all pairs of margins given \code{x[[icond]][ixcond]} are plotted where
    \code{icond} is taken from \code{x}.
  }
  \item{imargin}{vector of length 2.
    if given then conditional densities of the (\code{imargin[1]},
    \code{imargin[2]}) pair of margins given all values of \code{x[[icond]]}
    are plotted.
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
\seealso{
  \code{\link{NMixPredCondDensJoint2}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{dplot}
