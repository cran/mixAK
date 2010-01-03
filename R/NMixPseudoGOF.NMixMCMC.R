##
##  PURPOSE:   Pseudo goodness-of-fit test for a normal mixture
##             * method for class NMixMCMC
##
##  AUTHOR:    Arnost Komarek
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   20/08/2009
##
##  FUNCTIONS: NMixPseudoGOF.NMixMCMC.R
##
## ==================================================================

## *************************************************************
## NMixPseudoGOF.NMixMCMC
## *************************************************************
NMixPseudoGOF.NMixMCMC <- function(x, y, breaks, nbreaks=10, digits=3, ...)
{
  if (x$prior$priorK != "fixed") stop("only implemented for models with fixed number of components")

  return(NMixPseudoGOF.default(x=y, scale=x$scale, w=x$poster.mean.w, mu=x$poster.mean.mu, Sigma=x$poster.mean.Sigma, breaks=breaks, nbreaks=nbreaks, digits=digits))  
}
