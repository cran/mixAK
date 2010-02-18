##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2010
##
##  FUNCTIONS: tracePlots.NMixMCMC (16/01/2010)
##
## ======================================================================

## *************************************************************
## tracePlots.NMixMCMC
## *************************************************************
tracePlots.NMixMCMC <- function(x, param=c("Emix", "SDmix", "Cormix", "w", "mu", "sd", "gammaInv"),
                                relabel=FALSE, order,
                                auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", ...)
{
  param <- match.arg(param)

  ### Determine component in x where to look for required chains  
  if (param %in% c("Emix", "SDmix", "Cormix")) obj <- "mixture"
  else if (param == "sd") obj <- "Sigma"
       else               obj <- param

  if (param %in% c("w", "mu", "sd") & x$prior$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param %in% c("Emix", "SDmix")) nparam <- x$dim
  else if (param == "Cormix") nparam <- (x$dim * (x$dim + 1)) / 2 - x$dim
       else if (param == "w") nparam <- x$K[1]
            else if (param %in% c("mu", "sd")) nparam <- x$K[1] * x$dim
                 else if (param == "gammaInv") nparam <- x$dim
  if (!nparam){
    cat("\nNothing to plot.\n")
    return(invisible(x))    
  }  
    
  ### Plotting arguments
  if (length(col) == 1)  col <- rep(col, nparam)
  if (length(xlab) == 1) xlab <- rep(xlab, nparam)
  if (length(main) == 1) main <- rep(main, nparam)  
  if (length(col) != nparam)  stop("col must be of length", nparam)
  if (length(xlab) != nparam) stop("xlab must be of length", nparam)
  if (length(main) != nparam) stop("main must be of length", nparam)    

  if (!missing(ylab)){
    if (length(ylab) == 1)  ylab <- rep(ylab, nparam)
    if (length(ylab) != nparam) stop("ylab must be of length", nparam)    
  }
  
  ### Layout
  if (auto.layout){
    oldPar <- par(bty="n")
    layout(autolayout(nparam))
    on.exit(oldPar)
  }  

  ### Iteration index
  itIndex <- (x$nMCMC["burn"] + 1):(x$nMCMC["burn"] + x$nMCMC["keep"])

  ### Traceplots if related to moments of the mixture
  if (param %in% c("Emix", "SDmix", "Cormix")){
    if (param == "Emix") COLS <- paste("y.Mean.", 1:nparam, sep="")
    else if (param == "SDmix") COLS <- paste("y.SD.", 1:nparam, sep="")
         else if (param == "Cormix"){
           Imat <- diag(x$dim)
           rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
           colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)] 
           COLS <- paste("y.Corr.", rowsI, ".", colsI, sep="")           
         }
    
    if (missing(ylab)) ylab <- COLS

    for (i in 1:nparam){
      plot(itIndex, x[[obj]][, COLS[i]], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
    }  
  }

  else{

    ### Traceplots of mixture weights, means or standard deviations (possibly re-labeled)
    if (param %in% c("w", "mu", "sd")){
      if (!relabel) order <- matrix(rep(1:x$K[1], nrow(x[[obj]])), ncol=x$K[1], byrow=TRUE)
      if (relabel & missing(order)) order <- x$order
      if (nrow(order) != nrow(x[[obj]])) stop("order has incompatible number of rows")
      if (ncol(order) != x$K[1])         stop("order has incompatible number of columns")
      if (any(apply(order, 1, min) != 1))      stop("order must contain 1 in each row")
      if (any(apply(order, 1, max) != x$K[1])) stop("order must contain", x$K[1], "in each row")

      if (param == "w"){
        if (missing(ylab)) ylab <- colnames(x[[obj]])

        for (k in 1:x$K[1]){
          plot(itIndex, x[[obj]][cbind(1:nrow(x[[obj]]), order[,k])],
               type="l", xlab=xlab[k], ylab=ylab[k], col=col[k], main=main[k], ...)
        }
      }
      else{
        if (param == "mu"){     ### Draw traceplots of shifted and scaled mixture means
          if (missing(ylab)) ylab <- colnames(x[[obj]])

          for (k in 1:x$K[1]){
            for (j in 1:x$dim){
              i <- (k-1)*x$dim + j
              plot(itIndex, x$scale$shift[j] + x$scale$scale[j]*x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*x$dim + j)],
                   type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
            }
          }
        }
        else{
          if (param == "sd"){   ### Draw traceplots of scaled mixture standard deviations
            if (missing(ylab)) ylab <- paste("sd.", rep(1:x$K[1], each=x$dim), ".", rep(1:x$dim, x$K[1]), sep="")
            LTp <- (x$dim * (x$dim + 1)) / 2
            Idiag <- matrix(0, nrow=x$dim, ncol=x$dim)
            Idiag[lower.tri(Idiag, diag=TRUE)] <- 1:LTp
            jdiag <- diag(Idiag)
        
            for (k in 1:x$K[1]){
              for (j in 1:x$dim){
                i <- (k-1)*x$dim + j
                plot(itIndex, x$scale$scale[j]*sqrt(x[[obj]][cbind(1:nrow(x[[obj]]), (order[,k]-1)*LTp + jdiag[j])]),
                     type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
              }
            }        
          }  
        }  
      }  
    }

    ### Traceplots of all other parameters
    else{
      if (missing(ylab)) ylab <- colnames(x[[obj]])
      
      for (i in 1:nparam){
        plot(itIndex, x[[obj]][, i], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
      }
    }  
  }    
  
  return(invisible(x))
}  
