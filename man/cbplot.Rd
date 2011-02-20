\name{cbplot}
\alias{cbplot}
\title{
  Plot a function together with its confidence/credible bands
}
\description{
  This routine typically plots a function together with its confidence
  or credible bands. The credible band can be indicated either by
  additional lines or by a shaded region or by both.
}
\usage{
cbplot(x, y, low, upp, band.type=c("ls", "s", "l"), add=FALSE,
       col="darkblue", lty=1, lwd=2,
       cbcol=col, cblty=4, cblwd=lwd,
       scol=rainbow_hcl(1, start=180, end=180), slwd=5,
       xlim, ylim, xlab, ylab, main="", sub="", \dots)
}
\arguments{
  \item{x}{a numeric vector with \code{x} coordinates corresponding to
    \code{y}, \code{low}, \code{upp}.
  }
  \item{y}{a numeric vector with \code{y} coordinates of the function to
    plot.
  }
  \item{low}{a numeric vector with \code{y} coordinates of the lower
    limit of the credible band.
  }
  \item{upp}{a numeric vector with \code{y} coordinates of the upper
    limit of the credible band.
  }  
  \item{band.type}{a character which specifies the graphical way to show
    the credible band, \dQuote{ls} stands for line and shaded region,
    \dQuote{s} stands for shaded region only and \dQuote{l} stands for
    line only.
  }
  \item{add}{if \code{TRUE} then everything is added to the current
    plot.
  }
  \item{col, lty, lwd}{graphical paramters to draw the x-y line.
  }
  \item{cbcol, cblty, cblwd}{graphical parameters to draw the x-low and
    x-upp lines.
  }
  \item{scol, slwd}{graphical parameters for the shaded region between
    the credible/confidence bounds.
  }
  \item{xlim, ylim, xlab, ylab, main, sub}{other graphical parameters.
  }
  \item{\dots}{additional arguments passed to the \code{plot} function.
  }  
}
\value{
  \code{invisible(x)}
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### Artificial credible bands around the CDF's of N(100, 15*15)
### and N(80, 10*10)
iq <- seq(55, 145, length=100)
Fiq <- pnorm(iq, 100, 15)
low <- Fiq - 0.1
upp <- Fiq + 0.1

iq2 <- seq(35, 125, length=100)
Fiq2 <- pnorm(iq, 80, 10)
low2 <- Fiq2 - 0.1
upp2 <- Fiq2 + 0.1

cbplot(iq, Fiq, low, upp, xlim=c(35, 145))
cbplot(iq2, Fiq2, low2, upp2, add=TRUE, col="red4",
       scol=rainbow_hcl(1, start=20, end=20))
}
\keyword{dplot}

