##
##  PURPOSE:   Plot the function together with its confidence/credible band
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   08/05/2010
##      LOG:   02/06/2016 argument type added to allow drawing of stepwise functions
##                        (useful for Kaplan-Meier etc.)
##
##  FUNCTIONS: cbplot (08/05/2010)        
##             
## ======================================================================

## *************************************************************
## cbplot
## *************************************************************
##
cbplot <- function(x, y, low, upp, type=c("l", "s"), band.type=c("ls", "s", "l"), add=FALSE,
                   col="darkblue", lty=1, lwd=2,
                   cbcol=col, cblty=4, cblwd=lwd,
                   scol=rainbow_hcl(1, start=180, end=180), slwd=5,
                   xlim, ylim, xlab, ylab, main="", sub="", cex.lab=1, cex.axis=1, ...)
{
  type <- match.arg(type)
    
  if (length(x) < 2) stop("x must have the length of 2 or more")
  if (length(x) != length(y)) stop("incorrect y")
  if (length(x) != length(low)) stop("incorrect low")
  if (length(x) != length(upp)) stop("incorrect upp")

  if (missing(xlim)) xlim <- range(x, na.rm = TRUE)
  if (missing(ylim)) ylim <- range(c(y, low, upp), na.rm = TRUE)
  if (missing(xlab)) xlab <- substitute(x)
  if (missing(ylab)) ylab <- substitute(y)

  band.type <- match.arg(band.type)
  
  if (!add) plot(xlim, ylim, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub, cex.axis=cex.axis, cex.lab=cex.lab, ...)
  if (band.type %in% c("ls", "s")) segments(x0=x, y0=low, x1=x, y1=upp, col=scol, lwd=slwd, lend="round")
  lines(x, y, type=type, col=col, lwd=lwd, lty=lty)
  if (band.type %in% c("ls", "l")){
    lines(x, low, type=type, col=cbcol, lwd=cblwd, lty=cblty)
    lines(x, upp, type=type, col=cbcol, lwd=cblwd, lty=cblty)        
  }  
  
  return(invisible(x))  
}

