\name{plot.NMixPlugDensMarg}
\alias{plot.NMixPlugDensMarg}
\title{
  Plot computed marginal predictive densities
}
\description{
  This is a basic plotting tool to visualize computed marginal
  plug-in estimates of densities, see \code{\link{NMixPlugDensMarg}}.
}
\usage{
\method{plot}{NMixPlugDensMarg}(x, auto.layout=TRUE,
    type="l", col="darkblue", lty=1, lwd=1, main, xlab, ylab, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPlugDensMarg}.}
  \item{auto.layout}{if \code{TRUE} then the function determines itself
    how to divide the plotting region to draw densities for all margins.
  }
  \item{type}{type of the plot.}
  \item{col}{color used to draw the lines.}
  \item{lty}{type of the line.}
  \item{lwd}{line width.}
  \item{main}{main title of the plot. Either character which is
    replicated or a vector of characters of the length equal to the number of margins.}
  \item{xlab}{label for the x-axis. Either character which is
    replicated or a vector of characters of the length equal to the number of margins.}
  \item{ylab}{label for the y-axis. Either character which is
    replicated or a vector of characters of the length equal to the number of margins.}
  \item{\dots}{additional arguments passed to the \code{plot} function.} 
}
\value{
  \code{invisible(x)}
}
\seealso{
  \code{\link{NMixPlugDensMarg}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{dplot}
