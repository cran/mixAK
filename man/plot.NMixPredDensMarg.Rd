\name{plot.NMixPredDensMarg}
\alias{plot.NMixPredDensMarg}
\title{
  Plot computed marginal predictive densities
}
\description{
  This is a basic plotting tool to visualize computed marginal
  predictive densities, see \code{\link{NMixPredDensMarg}}.
}
\usage{
\method{plot}{NMixPredDensMarg}(x, K=0, auto.layout=TRUE,
    type="l", col="darkblue", lty=1, lwd=1, main, xlab, ylab, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPredDensMarg}.}
  \item{K}{if equal to \code{0} then the overall predictive densities
    are plotted taken from the \code{dens} part of the object \code{x}.

    If higher than \code{0} then the predictive density conditioned by
    the value of \code{K} is plotted (taken from the \code{densK} part
    of the object \code{x}).
  }
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
\references{
  Komárek, A. (2009).
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis}, \bold{53}, 3932-3947.
}
\seealso{
  \code{\link{NMixPredDensMarg}}, \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
## See additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Galaxy.pdf
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Faithful.pdf
## http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/Tandmob.pdf
##
}
\keyword{dplot}
