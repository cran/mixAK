\name{plot.NMixPredCondCDFMarg}
\alias{plot.NMixPredCondCDFMarg}
\title{
  Plot computed univariate conditional predictive cumulative
  distribution functions
}
\description{
  This is a basic plotting tool to visualize computed posterior
  predictive estimates
  of univariate conditional cdf's, see \code{\link{NMixPredCondCDFMarg}}.
}
\usage{
\method{plot}{NMixPredCondCDFMarg}(x, ixcond, imargin, prob, over=FALSE,
    auto.layout=TRUE, type="l", lwd=1, lty, col, qlwd=1, qlty, qcol,
    main, xlab, ylab, ylim,
    annot=TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPredCondCDFMarg}.}
  \item{ixcond}{if given then conditional cdf's of all margins given \code{x[[icond]][ixcond]} are plotted where
    \code{icond} is taken from \code{x}.
  }
  \item{imargin}{if given then conditional cdf's of the \code{imargin}-th margin given all values of \code{x[[icond]]}
    are plotted - either separately or all in one plot.
  }
  \item{prob}{probabilities of pointwise posterior quantiles which
    should be added to the plot. Computed values of requested posterior
    quantiles must be present in the object \code{x} (see argument
    \code{prob} of \code{\link{NMixPredCondCDFMarg}}).
  }
  \item{over}{logical. If \code{TRUE} and \code{imargin} is given then
    all conditional cdf's are drawn in one plot.
  }  
  \item{auto.layout}{if \code{TRUE} then the function determines itself
    how to divide the plotting region to draw the computed cdf's.
  }
  \item{type}{type of the plot.}
  \item{lwd}{line width.}  
  \item{lty}{type of the line. It can be a vector in which
    case different lines are drawn in different types.}
  \item{col}{color used to draw the lines. It can be a vector in which
    case different lines are drawn in different colors.}  
  \item{qlwd}{line width for pointwise posterior quantiles.}
  \item{qlty}{type of the line for pointwise posterior quantiles.}
  \item{qcol}{color used to draw pointwise posterior quantiles.}    
  \item{main}{main title of the plot. Either character which is
    replicated or a vector of characters.}
  \item{xlab}{label for the x-axis. Either character which is
    replicated or a vector of characters.}
  \item{ylab}{label for the y-axis. Either character which is
    replicated or a vector of characters.}
  \item{ylim}{limits for the y-axis.}
  \item{annot}{if \code{TRUE} and \code{imargin} is given and
    \code{over} is \code{TRUE} then a legend is added to the plot.}
  \item{\dots}{additional arguments passed to the \code{plot} function.} 
}
\value{
  \code{invisible(x)}
}
\seealso{
  \code{\link{NMixPredCondCDFMarg}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{dplot}
