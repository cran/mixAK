##
##  PURPOSE:   Convert fitted distribution of Y=log(T) into distribution of T
##             * method for objects of class NMixPredCDFMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: Y2T.NMixPredCDFMarg (09/06/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredCDFMarg
## *************************************************************
Y2T.NMixPredCDFMarg <- function(x, ...)
{
  ### Unlog the grids and get T = exp(Y)
  x$x <- lapply(x$x, exp)

  return(x)
}  
