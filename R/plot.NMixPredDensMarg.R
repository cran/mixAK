##
##  PURPOSE:   Plotting of computed predictive marginal (univariate) densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##
##  FUNCTIONS: plot.NMixPredDensMarg (03/12/2007)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredDensMarg
## *************************************************************
plot.NMixPredDensMarg <- function(x, K=0, auto.layout=TRUE, type="l", col="darkblue", lty=1, lwd=1, main, xlab, ylab, ...)
{
  p <- length(x$x)  
  nfig <- p

  miss.main <- missing(main)
  
  if (missing(xlab)){
    if (p == 1) xlab <- "x"
    else        xlab <- paste("x", 1:p, sep="")
  }else{
    if (length(xlab) == 1){
      xlab <- rep(xlab, p)
    }else{
      if (length(xlab) != p) stop(paste("xlab must be of length 1 or ", p, sep=""))
    }  
  }  
  
  if (missing(ylab)){
    ylab <- rep("Density", p)
  }else{
    if (length(ylab) == 1){
      ylab <- rep(ylab, p)
    }else{
      if (length(ylab) != p) stop(paste("ylab must be of length 1 or ", p, sep=""))
    }  
  }  

  #if (length(K) != 1) stop("K must be of length 1")
  if (any(K < 0)) stop("K must not be negative")
  if (any(K > length(x$densK[[1]]))) stop("K is too high")

  percK <- paste(" (", round(x$propK*100, 1), "%)", sep="")
      
  if (auto.layout){    
    LAY <- autolayout(p)    
    oldPar <- par(bty="n")
    layout(LAY)
    on.exit(oldPar)
  }

  if (length(K) == 1){
    for (m0 in 1:p){
      if (K == 0){
        dx <- x$dens[[paste(m0)]]
        main2  <- paste(" (MCMC length = ", x$MCMC.length, ")", sep="")
        main2b <- paste("MCMC length = ", x$MCMC.length, sep="")                  
      }else{
        dx <- x$densK[[paste(m0)]][[K]]
        main2  <- paste(",  K = ", K, percK[K], sep="")
        main2b <- paste("K = ", K, percK[K], sep="")                  
      }
      
      if (miss.main){
        if (p == 1) MAIN <- main2b
        else        MAIN <- paste("Margin ", m0, main2, sep="")
      }else{
        if (length(main) == 1) MAIN <- main
        else{
          if (length(main) != p) stop(paste("main must be of length 1 or ", p, sep=""))
          MAIN <- main[m0]
        }  
      }  

      plot(x$x[[m0]], dx, type=type, col=col[1], main=MAIN, xlab=xlab[m0], ylab=ylab[m0], lty=lty[1], lwd=lwd[1], ...)
    }    
  }else{
    if (length(col) == 1) col <- rep(col, length(K))
    if (length(col) != length(K)) stop("incorrect length(col) argument")
    if (length(lty) == 1) lty <- rep(lty, length(K))
    if (length(lty) != length(K)) stop("incorrect length(lty) argument")
    if (length(lwd) == 1) lwd <- rep(lwd, length(K))
    if (length(lwd) != length(K)) stop("incorrect length(lwd) argument")        
    for (m0 in 1:p){
      if (miss.main){
        if (p == 1) main <- paste("MCMC length = ", x$MCMC.length, sep="")
        else        main <- paste("Margin ", m0, " (MCMC length = ", x$MCMC.length, ")", sep="")
      }else{
        if (length(main) == 1) MAIN <- main
        else{
          if (length(main) != p) stop(paste("main must be of length 1 or ", p, sep=""))
          MAIN <- main[m0]
        }
      }  
      
      for (kk in 1:length(K)){
        if (K[kk] == 0) dx <- x$dens[[paste(m0)]]
        else            dx <- x$densK[[paste(m0)]][[K[kk]]]
        if (kk == 1) plot(x$x[[m0]], dx, type=type, col=col[1], main=main, xlab=xlab[m0], ylab=ylab[m0], lty=lty[1], lwd=lwd[1], ...)
        else         lines(x$x[[m0]], dx, col=col[kk], lty=lty[kk], lwd=lwd[kk])
      } 
    }            
  }  
    
  return(invisible(x))   
}  


