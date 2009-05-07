##
##  PURPOSE:   Plotting of computed predictive pairwise joint densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##
##  FUNCTION:  plot.NMixPredDensJoint2 (03/12/2007)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredDensJoint2
## *************************************************************
plot.NMixPredDensJoint2 <- function(x, K=0, auto.layout=TRUE, col="darkblue", lwd=1, main, xlab, ylab, ...)
{
  p <- length(x$x)  
  nfig <- p * (p-1)/2

  miss.main <- missing(main)
  miss.xlab <- missing(xlab)
  miss.ylab <- missing(ylab)
  
  if (length(K) != 1) stop("K must be of length 1")
  if (K < 0) stop("K must not be negative")
  if (K > length(x$densK[[1]])) stop("K is too high")

  percK <- paste(" (", round(x$propK*100, 1), "%)", sep="")
  
  if (p == 2){
    if (auto.layout){
      oldPar <- par(mfrow=c(1, 1), bty="n")
      on.exit(oldPar)
    }  

    if (K == 0){
      dx <- x$dens[["1-2"]]
      main2 <- paste(" (MCMC length = ", x$MCMC.length, ")", sep="")
    }else{
      dx <- x$densK[["1-2"]][[K]]
      main2 <- paste(",  K = ", K, percK[K], sep="")
    }  

    if (miss.main) main <- paste("Margins (", 1, ", ", 2, ")", main2, sep="")
    if (miss.xlab) xlab <- "x1"
    if (miss.ylab) ylab <- "x2"
    contour(x$x[[1]], x$x[[2]], dx, col=col, main=main, xlab=xlab, ylab=ylab, lwd=lwd, ...)    
    
  }else{
    if (auto.layout){
      if (p == 3) lay <- c(2, 2)                                           ## 3 figures
      else if (p == 4) lay <- c(2, 3)                                      ## 6 figures
           else if (p == 5) lay <- c(3, 4)                                 ## 10 figures
                else if (p == 6) lay <- c(3, 5)                            ## 15 figures
                     else if (p == 7) lay <- c(4, 6)                       ## 21 figures
                          else if (p == 8) lay <- c(4, 7)                  ## 28 figures
                               else stop("layout must be given for p > 8")
      
      oldPar <- par(mfrow=lay, bty="n")
      on.exit(oldPar)
    }  
    
    for (m0 in 1:(p-1)){
      for (m1 in (m0+1):p){
        if (K == 0){
          dx <- x$dens[[paste(m0, "-", m1, sep="")]]
          main2 <- paste(" (MCMC length = ", x$MCMC.length, ")", sep="")          
        }else{
          dx <- x$densK[[paste(m0, "-", m1, sep="")]][[K]]
          main2 <- paste(",  K = ", K, percK[K], sep="")          
        }

        if (miss.main) main <- paste("Margins (", m0, ", ", m1, ")", main2, sep="")
        if (miss.xlab) xlab <- paste("x", m0, sep="")
        if (miss.ylab) ylab <- paste("x", m1, sep="")
        contour(x$x[[m0]], x$x[[m1]], dx, col=col, main=main, xlab=xlab, ylab=ylab, lwd=lwd, ...)
      }  
    }    
  }  
       
  return(invisible(x))   
}  


