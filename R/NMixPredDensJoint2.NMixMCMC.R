##
##  PURPOSE:   Computation of the predictive pairwise joint densities
##             * method for NMixMCMC objects
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##
##  FUNCTIONS: NMixPredDensJoint2.NMixMCMC (02/01/2008)
##             
## ======================================================================

## *************************************************************
## NMixPredDensJoint2.NMixMCMC
## *************************************************************
NMixPredDensJoint2.NMixMCMC <- function(x, grid, lgrid=50, scaled=FALSE, ...)
{
  if (missing(grid)){
    grid <- list()
    if (scaled){
      if (x$dim == 1){
        rangeGrid <- x$summ.z.Mean["Median"] + c(-3.5, 3.5)*x$summ.z.SDCorr["Median"]
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dim){
          rangeGrid <- x$summ.z.Mean["Median", i] + c(-3.5, 3.5)*x$summ.Z.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }else{
      if (x$dim == 1){
        rangeGrid <- x$summ.y.Mean["Median"] + c(-3.5, 3.5)*x$summ.y.SDCorr["Median"]
        grid[[1]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)      
      }else{     
        for (i in 1:x$dim){
          rangeGrid <- x$summ.y.Mean["Median", i] + c(-3.5, 3.5)*x$summ.y.SDCorr["Median", (i-1)*(2*x$dim - i + 2)/2 + 1]
          grid[[i]] <- seq(rangeGrid[1], rangeGrid[2], length=lgrid)
        }
      }  
    }
    names(grid) <- paste("x", 1:x$dim, sep="")    
  }
  
  if (x$dim == 1) if (is.numeric(grid)) grid <- list(x1=grid)
  if (!is.list(grid)) stop("grid must be a list")

  if (scaled) scale <- list(shift=0, scale=1)  
  else        scale <- x$scale

  if (x$prior$priorK == "fixed"){
    return(NMixPredDensJoint2.default(x=grid, scale=scale, K=x$K, w=as.numeric(t(x$w)), mu=as.numeric(t(x$mu)), Li=as.numeric(t(x$Li)), Krandom=FALSE))
  }else{
    return(NMixPredDensJoint2.default(x=grid, scale=scale, K=x$K, w=x$w, mu=x$mu, Li=x$Li, Krandom=TRUE))
  }    
}  




