##
##  PURPOSE:   Reversible jump MCMC for a normal mixture model
##             * print method for the resulting object
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   29/10/2007
##
##  FUNCTION:  print.NMixMCMC 
##             
##
## ======================================================================

## *************************************************************
## print.NMixMCMC
## *************************************************************
print.NMixMCMC <- function(x, ...)
{
  cat("\n")
  if (x$prior$priorK == "fixed"){
    cat(paste("     ", x$prior$Kmax, " component normal mixture estimated using MCMC\n", sep=""))
    cat("     ================================================\n")
  }else{
    cat(paste("     Normal mixture with at most ", x$prior$Kmax, " components estimated using RJ-MCMC\n", sep=""))
    cat("     =================================================================")
    cat("\nPosterior distribution of K:")
    cat("\n----------------------------")
    print(x$propK, ...)    
  }  

  cat("\nDeviance information criteria:")
  cat("\n------------------------------\n")
  print(x$DIC, ...)    
  
  cat("\nPosterior summary statistics for moments of mixture for original data:")
  cat("\n----------------------------------------------------------------------")
  if (x$dim == 1){
    cat("\nMean:\n")
    print(x$summ.y.Mean, ...)

    cat("\nStandard deviation:\n")
    print(x$summ.y.SDCorr, ...)
  }else{
    cat("\nMeans:\n")
    print(x$summ.y.Mean, ...)

    cat("\nStandard deviations and correlations:\n")
    print(x$summ.y.SDCorr, ...)
  }  

  if (FALSE) {
    cat("\nPosterior summary statistics for moments of mixture for scaled data:")
    cat("\n--------------------------------------------------------------------")
    if (x$dim == 1){
      cat("\nMean:\n")
      print(x$summ.z.Mean, ...)

      cat("\nStandard deviation:\n")
      print(x$summ.z.SDCorr, ...)
    }else{
      cat("\nMeans:\n")
      print(x$summ.z.Mean, ...)

      cat("\nStandard deviations and correlations:\n")
      print(x$summ.z.SDCorr, ...)
    }
  }  

  return(invisible(x))
}
