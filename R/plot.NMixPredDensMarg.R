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
  miss.xlab <- missing(xlab)
  miss.ylab <- missing(ylab)
  
  #if (length(K) != 1) stop("K must be of length 1")
  if (any(K < 0)) stop("K must not be negative")
  if (any(K > length(x$densK[[1]]))) stop("K is too high")

  percK <- paste(" (", round(x$propK*100, 1), "%)", sep="")
      
  if (auto.layout){
    if (p == 1) lay <- c(1, 1)
    else if (p <= 2) lay <- c(1, 2)
         else if (p <= 4) lay <- c(2, 2)
              else if (p <= 6) lay <- c(2, 3)
                   else if (p <= 9) lay <- c(3, 3)
                        else if (p <= 12) lay <- c(3, 4)
                             else if (p <= 16) lay <- c(4, 4)
                                  else if (p <= 20) lay <- c(4, 5)
                                  else stop("layout must be given for p > 20")
    
    oldPar <- par(mfrow=lay, bty="n")
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
        if (p == 1) main <- main2b
        else        main <- paste("Margin ", m0, main2, sep="")
      }  
      if (miss.xlab) if (p == 1) xlab <- "x" else xlab <- paste("x", m0, sep="")
      if (miss.ylab) ylab <- "Density"
      plot(x$x[[m0]], dx, type=type, col=col[1], main=main, xlab=xlab, ylab=ylab, lty=lty[1], lwd=lwd[1], ...)
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
      }
      if (miss.xlab) if (p == 1) xlab <- "x" else xlab <- paste("x", m0, sep="")
      if (miss.ylab) ylab <- "Density"      
      for (kk in 1:length(K)){
        if (K[kk] == 0) dx <- x$dens[[paste(m0)]]
        else            dx <- x$densK[[paste(m0)]][[K[kk]]]
        if (kk == 1) plot(x$x[[m0]], dx, type=type, col=col[1], main=main, xlab=xlab, ylab=ylab, lty=lty[1], lwd=lwd[1], ...)
        else         lines(x$x[[m0]], dx, col=col[kk], lty=lty[kk], lwd=lwd[kk])
      } 
    }            
  }  
    
  return(invisible(x))   
}  


