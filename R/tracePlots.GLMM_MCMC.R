##
##  PURPOSE:   Drawing traceplots of selected parameters from fitted objects
##             * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2010
##
##  FUNCTIONS: tracePlots.GLMM_MCMC (16/01/2010)
##
## ======================================================================

## *************************************************************
## tracePlots.GLMM_MCMC
## *************************************************************
tracePlots.GLMM_MCMC <- function(x, param=c("beta", "Eb", "SDb", "Corb", "sigma_eps", "w_b", "mu_b", "gammaInv_b", "gammaInv_eps"),
                                 auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", ...)
{
  param <- match.arg(param)
  obj <- ifelse(param %in% c("Eb", "SDb", "Corb"), "mixture_b", param)      ### component in x where to look for required chains

  if (param %in% c("w_b", "mu_b") & x$prior.b$priorK != "fixed") stop("Not implemented for this value of param.")
  
  ### Number of parameters to plot
  if (param == "beta") nparam <- sum(x$p)
  else if (param %in% c("Eb", "SDb")) nparam <- x$dimb
       else if (param == "Corb") nparam <- (x$dimb * (x$dimb + 1)) / 2 - x$dimb
            else if (param == "sigma_eps") nparam <- x$R["Rc"]
                 else if (param == "w_b") nparam <- x$K_b[1]
                      else if (param == "mu_b") nparam <- x$K_b[1] * x$dimb
                           else if (param == "gammaInv_b") nparam <- x$dimb
                                else if (param == "gammaInv_eps") nparam <- x$R["Rc"]
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
  if (param %in% c("Eb", "SDb", "Corb")){
    if (param == "Eb") COLS <- paste("b.Mean.", 1:nparam, sep="")
    else if (param == "SDb") COLS <- paste("b.SD.", 1:nparam, sep="")
         else if (param == "Corb"){
           Imat <- diag(x$dimb)
           rowsI <- row(Imat)[lower.tri(row(Imat), diag=FALSE)]
           colsI <- col(Imat)[lower.tri(col(Imat), diag=FALSE)] 
           COLS <- paste("b.Corr.", rowsI, ".", colsI, sep="")           
         }
    
    if (missing(ylab)) ylab <- COLS
    
    for (i in 1:nparam){
      plot(itIndex, x[[obj]][, COLS[i]], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
    }  
  }

  ### Traceplots of all other parameters
  else{
    if (missing(ylab)) ylab <- colnames(x[[obj]])

    if (param == "mu_b"){     ### Draw traceplots of shifted and scaled mixture means
      for (k in 1:x$K_b[1]){
        for (j in 1:x$dimb){
          i <- (k-1)*x$dimb[1] + j
          plot(itIndex, x$scale.b$shift[j] + x$scale.b$scale[j]*x[[obj]][, paste("mu.", k, ".", j, sep="")],
               type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
        }
      }
    }
    else{    
      for (i in 1:nparam){
        plot(itIndex, x[[obj]][, i], type="l", xlab=xlab[i], ylab=ylab[i], col=col[i], main=main[i], ...)
      }
    }  
  }  
  
  return(invisible(x))
}


