\name{plot.NMixPredCondDensMarg}
\alias{plot.NMixPredCondDensMarg}
\title{
  Plot computed univariate conditional predictive densities
}
\description{
  This is a basic plotting tool to visualize computed plug-in estimates
  of univariate conditional densities, see \code{\link{NMixPredCondDensMarg}}.
}
\usage{
\method{plot}{NMixPredCondDensMarg}(x, ixcond, imargin, over=FALSE,
    auto.layout=TRUE, type="l", lwd=1, lty, col, main, xlab, ylab, ylim,
    annot=TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixPredCondDensMarg}.}
  \item{ixcond}{if given then conditional densities of all margins given \code{x[[icond]][ixcond]} are plotted where
    \code{icond} is taken from \code{x}.
  }
  \item{imargin}{if given then conditional densities of the \code{imargin}-th margin given all values of \code{x[[icond]]}
    are plotted - either separately or all in one plot.
  }
  \item{over}{logical. If \code{TRUE} and \code{imargin} is given then
    all conditional densities are drawn in one plot.
  }  
  \item{auto.layout}{if \code{TRUE} then the function determines itself
    how to divide the plotting region to draw the computed densities.
  }
  \item{type}{type of the plot.}
  \item{lwd}{line width.}  
  \item{col}{color used to draw the lines. It can be a vector in which
    case different lines are drawn in different colors.}
  \item{lty}{type of the line. It can be a vector in which
    case different lines are drawn in different types.}
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
  \code{\link{NMixPredCondDensMarg}}, \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{dplot}
