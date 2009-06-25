##
##  PURPOSE:   Plotting of computed pairwise bivariate conditional densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   29/05/2009
##
##  FUNCTIONS: plot.NMixPlugCondDensJoint2 (29/05/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPlugCondDensJoint2
## *************************************************************
##
plot.NMixPlugCondDensJoint2 <- function(x, ixcond, imargin, contour=FALSE, auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  return(plot.NMixPredCondDensJoint2(x=x, ixcond=ixcond, imargin=imargin, contour=contour, auto.layout=auto.layout,
                                     col=col, lwd=lwd, main=main, xylab=xylab, ...))

}  
