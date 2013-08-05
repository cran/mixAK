##
##  PURPOSE:   Re-labeling of the MCMC output.
##             * method for objects of class GLMM_MCMClist
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   25/06/2013
##
##  FUNCTION:  NMixRelabel.GLMM_MCMClist (25/06/2013) 
##
## ======================================================================

## *************************************************************
## NMixRelabel.GLMM_MCMClist
## *************************************************************
NMixRelabel.GLMM_MCMClist <- function(object, type = c("mean", "weight", "stephens"), par,
                                      prob = c(0.025, 0.5, 0.975), keep.comp.prob = FALSE, info, silent = FALSE,
                                      parallel = FALSE, ...)
{
  if (parallel){
    RAlg <- NMixRelabelAlgorithm(type = type, par = par, dim = object[[1]]$dimb)
    if (missing(info)) info <- object[[1]]$nMCMC["keep"]
    
    #require("parallel")
    if (parallel::detectCores() < 2) warning("It does not seem that at least 2 CPU cores are available needed for efficient parallel re-labelling of the two chains.")      
    cl <- parallel::makeCluster(2)      
    if (!silent){
      cat("\nParallel re-labelling of the two chains\n")
      cat("=======================================\n\n")
    }      
    tmpObj <- list(object[[1]], object[[2]])
    class(tmpObj) <- "GLMM_MCMClist"
    tmpObj <- parallel::parLapply(cl, tmpObj, NMixRelabel.GLMM_MCMC,
                                  type = RAlg$relabel$type, par = RAlg$relabel$par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)        
    parallel::stopCluster(cl)

    elemObj <- names(object)[-(1:2)]
    for (i in 1:length(elemObj)) tmpObj[[elemObj[i]]] <- object[[elemObj[i]]]
    
    return(tmpObj)
  }else{
    if (!silent){
      cat("\nRe-labelling chain number 1\n")
      cat("===========================\n")
    }  
    object[[1]] <- NMixRelabel(object[[1]], type = type, par = par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)

    if (!silent){    
      cat("\nRe-labelling chain number 2\n")
      cat("===========================\n\n")
    }  
    object[[2]] <- NMixRelabel(object[[2]], type = type, par = par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)
    cat("\n")
    
    return(object)
  }
}  
