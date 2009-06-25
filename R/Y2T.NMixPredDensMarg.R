##
##  PURPOSE:   Convert fitted distribution of Y=log(T) into distribution of T
##             * method for objects of class NMixPredDensMarg
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS: Y2T.NMixPredDensMarg (09/06/2009)
##
## ======================================================================

## *************************************************************
## Y2T.NMixPredDensMarg
## *************************************************************
Y2T.NMixPredDensMarg <- function(x, ...)
{
  ### Unlog the grids and get T = exp(Y)
  x$x <- lapply(x$x, exp)

  ### Multiply each marginal density by appropriate jacobian (t^{-1})
  for (i in 1:length(x$dens)){
    x$dens[[i]] <- x$dens[[i]] / x$x[[i]]
    for (k in 1:length(x$densK[[i]])){
      x$densK[[i]][[k]] <- x$densK[[i]][[k]] / x$x[[i]]
    }  
  }

  return(x)
}  
