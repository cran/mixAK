##
##  PURPOSE:   Plotting of computed predictive pairwise joint densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##  LOG:       25/05/2009:  xlab, ylab arguments removed
##                          xylab argument added instead
##                          contour argument added, function draws images by default
##             08/11/2011:  add.contour, col.add.contour arguments added
##
##  FUNCTION:  plot.NMixPredDensJoint2 (03/12/2007)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredDensJoint2
## *************************************************************
plot.NMixPredDensJoint2 <- function(x, K=0, contour=FALSE, add.contour=TRUE, col.add.contour="brown", auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  p <- length(x$x)  
  nfig <- p * (p-1)/2
  
  miss.main <- missing(main)
  miss.xylab <- missing(xylab)
  if (!miss.xylab){
    if (length(xylab) != p) stop("xylab must be of length", p)
  }  
  
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
    if (miss.xylab){
      xlab <- "x1"
      ylab <- "x2"
    }else{
      xlab <- xylab[1]
      ylab <- xylab[2]
    }
    if (contour){
      if (missing(col)) col <- "darkblue"      
      contour(x$x[[1]], x$x[[2]], dx, col=col, main=main, xlab=xlab, ylab=ylab, lwd=lwd, ...)
    }else{
      if (missing(col)){
        #require("colorspace")
        col <- rev(heat_hcl(33, c.=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))
      }
      image(x$x[[1]], x$x[[2]], dx, col=col, main=main, xlab=xlab, ylab=ylab, ...)
      if (add.contour) contour(x$x[[1]], x$x[[2]], dx, col=col.add.contour, lwd=lwd, add=TRUE)      
    }      
  }else{
    if (auto.layout){
      if (p == 3) LAY <- autolayout(3)                                          ## 3 figures
      else if (p == 4) LAY <- autolayout(6)                                     ## 6 figures
           else if (p == 5) LAY <- autolayout(10)                               ## 10 figures
                else if (p == 6) LAY <- autolayout(15)                          ## 15 figures
                     else if (p == 7) LAY <- autolayout(21)                     ## 21 figures
                          else if (p == 8) LAY <- autolayout(28)                ## 28 figures
                               else stop("layout must be given for p > 8")
      
      oldPar <- par(bty="n")
      layout(LAY)
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
        if (miss.xylab){
          xlab <- paste("x", m0, sep="")
          ylab <- paste("x", m1, sep="")
        }else{
          xlab <- xylab[m0]
          ylab <- xylab[m1]
        }

        if (contour){
          if (missing(col)) col <- "darkblue"          
          contour(x$x[[m0]], x$x[[m1]], dx, col=col, main=main, xlab=xlab, ylab=ylab, lwd=lwd, ...)
        }else{
          if (missing(col)){
            #require("colorspace")
            col <- rev(heat_hcl(33, c.=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))
          }
          image(x$x[[m0]], x$x[[m1]], dx, col=col, main=main, xlab=xlab, ylab=ylab, ...)
          if (add.contour) contour(x$x[[m0]], x$x[[m1]], dx, col=col.add.contour, lwd=lwd, add=TRUE)
        }  
      }  
    }    
  }  
       
  return(invisible(x))   
}  


