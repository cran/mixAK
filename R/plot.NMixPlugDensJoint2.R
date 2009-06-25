##
##  PURPOSE:   Plotting of computed pairwise joint densities (plug-in version)
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   28/05/2009
##
##  FUNCTION:  plot.NMixPlugDensJoint2 (28/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugDensJoint2
## *************************************************************
plot.NMixPlugDensJoint2 <- function(x, contour=FALSE, auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  return(plot.NMixPredDensJoint2(x=x, K=0, contour=contour, auto.layout=auto.layout, col=col, lwd=lwd, main=main, xylab=xylab, ...))
}
