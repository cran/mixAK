##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of continuous response and random effects
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    06/07/2009
##              03/08/2009:  version for continuous responses working
##
##  FUNCTIONS:  GLMM_MCMC
##
## ================================================================================================

## *************************************************************
## GLMM_MCMC
## *************************************************************
##
GLMM_MCMC <- function(y, dist="gaussian", id, x, z, random.intercept,
                      prior.beta, init.beta,                      
                      scale.b,    prior.b,   init.b,
                      prior.eps,  init.eps,
                      nMCMC=c(burn=10, keep=10, thin=1, info=10),
                      store=c(b=FALSE), keep.chains=TRUE)
{
  require("lme4")
  thispackage <- "mixAK"

  DEBUG <- FALSE
  
  parallel <- FALSE     ### at this moment, parallel computation not fully implemented

########## ========== Data ========== ##########
########## ========================== ##########
  dd <- GLMM_MCMCdata(y=y, dist=dist, id=id, x=x, z=z, random.intercept=random.intercept)
  rm(list=c("y", "dist", "id", "x", "z", "random.intercept"))
     ### use dd$y, dd$dist, dd$id, dd$x, dd$z, dd$random.intercept instead
     ### REMARK:  dd$x, dd$z are still without intercept column
  
  
########## ========== Initial fits ======================================== ##########
########## ========== and design information to be passed to C++ ========== ##########
########## ================================================================ ##########
  ifit <- GLMM_MCMCifit(do.init=TRUE, na.complete=FALSE,
                        y=dd$y, dist=dd$dist, id=dd$id, time=dd$time, x=dd$x, z=dd$z, random.intercept=dd$random.intercept,
                        xempty=dd$xempty, zempty=dd$zempty, Rc=dd$Rc, Rd=dd$Rd,
                        p=dd$p, p_fi=dd$p_fi, q=dd$q, q_ri=dd$q_ri, lbeta=dd$lbeta, dimb=dd$dimb)
  dd$x <- NULL
  dd$z <- NULL
     ### use ifit$x, ifit$z instead
     ### REMARK:  ifit$x, ifit$z contain intercept columns as well
  

########## ========== Prior distribution for fixed effects (beta) ========== ##########
########## ================================================================= ##########
  pbeta <- GLMM_MCMCprior.beta(prior.beta=prior.beta, lbeta=dd$lbeta)
  prior.beta <- pbeta$prior.beta
  pbeta$prior.beta <- NULL


########## ========== Prior distribution for error terms of gaussian responses ========== ##########
########## ============================================================================== ##########
  peps <- GLMM_MCMCprior.eps(prior.eps=prior.eps, Rc=dd$Rc, isigma=ifit$isigma, is.sigma=ifit$is.sigma)
  prior.eps <- peps$prior.eps
  peps$prior.eps <- NULL

  
########## ========== Shift and scale for random effects ========== ##########
########## ======================================================== ##########
  scb <- GLMM_MCMCscale.b(scale.b=scale.b, dimb=dd$dimb, iEranefVec=ifit$iEranefVec, iSDranefVec=ifit$iSDranefVec)
  scale.b <- scb$scale.b
  scb$scale.b <- NULL
  
  
########## ========== Prior distribution for random effects ========== ##########
########## =========================================================== ##########
  pbb <- GLMM_MCMCprior.b(prior.b=prior.b, scale.b=scale.b, dimb=dd$dimb, iEranefVec=ifit$iEranefVec, iSDranefVec=ifit$iSDranefVec)
  prior.b <- pbb$prior.b
  pbb$prior.b <- NULL
     
  
########## ========== Initial values for fixed effects (beta) ============== ##########
########## ================================================================= ##########  
  if (dd$lbeta){
    if (missing(init.beta)) init.beta <- ifit$ibeta
    if (!is.numeric(init.beta)) stop("init.beta must be numeric")
    if (length(init.beta) == 1) init.beta <- rep(init.beta, dd$lbeta)
    if (length(init.beta) != dd$lbeta) stop(paste("init.beta must be of length", dd$lbeta))
    if (is.null(names(init.beta))) names(init.beta) <- paste("beta", 1:dd$lbeta, sep="")
  }else{
    init.beta <- 0
  }  
  
  
########## ========== Initial values for parameters related to the distribution of error terms of gaussian responses ========== ##########
########## ==================================================================================================================== ##########
  if (dd$Rc){
    if (missing(init.eps)) init.eps <- list()
    if (!is.list(init.eps)) stop("init.eps must be a list")    
    
    ininit.eps <- names(init.eps)
    ieps.sigma    <- match("sigma", ininit.eps, nomatch=NA)
    ieps.gammaInv <- match("gammaInv", ininit.eps, nomatch=NA)    

    ##### init.eps:  sigma
    ##### -----------------------------------------------
    if (is.na(ieps.sigma)) init.eps$sigma <- ifit$isigma[ifit$is.sigma]
    if (length(init.eps$sigma) == 1) init.eps$sigma <- rep(init.eps$sigma, dd$Rc)
    if (length(init.eps$sigma) != dd$Rc) stop(paste("init.eps$sigma must be of length", dd$Rc))
    if (any(is.na(init.eps$sigma))) stop("NA in init.eps$sigma")        
    if (any(init.eps$sigma <= 0)) stop(paste("init.eps$sigma must be higher than ", 0, sep=""))
    names(init.eps$sigma) <- paste("sigma", 1:dd$Rc, sep="")

    ##### init.eps:  gammaInv
    ##### -----------------------------------------------
    if (is.na(ieps.gammaInv)) init.eps$gammaInv <- prior.eps$zeta * ifit$isigma[ifit$is.sigma]^2
    if (length(init.eps$gammaInv) == 1) init.eps$gammaInv <- rep(init.eps$gammaInv, dd$Rc)
    if (length(init.eps$gammaInv) != dd$Rc) stop(paste("init.eps$gammaInv must be of length", dd$Rc))
    if (any(is.na(init.eps$gammaInv))) stop("NA in init.eps$gammaInv")        
    if (any(init.eps$gammaInv <= 0)) stop(paste("init.eps$gammaInv must be higher than ", 0, sep=""))
    names(init.eps$gammaInv) <- paste("gammaInv", 1:dd$Rc, sep="")        
  }else{
    init.eps <- list(sigma=0, gammaInv=0)
  }  
  
  
########## ========== Initial values for parameters related to the distribution of random effects ========== ##########
########## ================================================================================================= ##########  
  if (dd$dimb){
    if (missing(init.b)) init.b <- list()
    if (!is.list(init.b)) stop("init.b must be a list")

    ininit.b    <- names(init.b)
    ib.b        <- match("b", ininit.b, nomatch=NA)
    ib.K        <- match("K", ininit.b, nomatch=NA)
    ib.w        <- match("w", ininit.b, nomatch=NA)
    ib.mu       <- match("mu", ininit.b, nomatch=NA)
    ib.Sigma    <- match("Sigma", ininit.b, nomatch=NA)
    ib.Li       <- match("Li", ininit.b, nomatch=NA)
    ib.gammaInv <- match("gammaInv", ininit.b, nomatch=NA)
    ib.r        <- match("r", ininit.b, nomatch=NA)

    ##### init.b:  b (not scaled and not shifted!!!)
    ##### ----------------------------------------------------  
    if (is.na(ib.b)) init.b$b <- as.matrix(ifit$ibMat)
    if (dd$dimb == 1) init.b$b <- matrix(as.numeric(init.b$b), ncol=1)
    if (is.data.frame(init.b$b)) init.b$b <- as.matrix(init.b$b)
    if (!is.matrix(init.b$b)) stop("init.b$b must be a matrix")
    if (ncol(init.b$b) != dd$dimb) stop(paste("init.b$b must have ", dd$dimb, " columns", sep=""))
    if (nrow(init.b$b) != ifit$I) stop(paste("init.b$b must have ", ifit$I, " rows", sep=""))
    if (is.null(rownames(init.b$b))) rownames(init.b$b) <- unique(dd$id)
    if (any(is.na(init.b$b))) stop("NA in init.b$b")
    
    ##### init.b:  K
    ##### ----------------------------------------------------  
    if (is.na(ib.K)){
      if (prior.b$priorK == "fixed") init.b$K <- prior.b$Kmax
      else                           init.b$K <- 1
    }
    if (prior.b$priorK == "fixed") init.b$K <- prior.b$Kmax
    if (length(init.b$K) != 1) stop("init.b$K must be of length 1")
    if (is.na(init.b$K)) stop("NA in init.b$K")
    if (init.b$K <= 0 | init.b$K > prior.b$Kmax) stop("init.b$K out of the range")
    
    ##### init.b:  w
    ##### ----------------------------------------------------  
    if (is.na(ib.w)){
      init.b$w <- rep(1, init.b$K)/init.b$K
    }  
    init.b$w <- as.numeric(init.b$w)
    if (length(init.b$w) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$w <- init.b$w[1:init.b$K]  
    names(init.b$w) <- paste("w", 1:init.b$K, sep="")
    if (any(is.na(init.b$w))) stop("NA in init.b$w")  
    if (length(init.b$w) != init.b$K) stop(paste("init.b$w must be of length ", init.b$K, sep=""))
    if (any(init.b$w < 0)) stop("init.b$w may not be negative")
    init.b$w <- init.b$w / sum(init.b$w)

    ##### init.b:  mu
    ##### ----------------------------------------------------
    Rbb <- 6 * (ifit$iSDranefVec / scale.b$scale)
    bbmin <- (ifit$iEranefVec - scale.b$shift) / scale.b$scale - 0.5 * Rbb
    bbmax <- (ifit$iEranefVec - scale.b$shift) / scale.b$scale + 0.5 * Rbb
    
    if (is.na(ib.mu)){      
      if (dd$dimb == 1){
        afstand <- Rbb/(init.b$K + 1)
        init.b$mu <- seq(bbmin+afstand, bbmax-afstand, length=init.b$K)
      }else{
        afstand <- Rbb/(init.b$K + 1)
        init.b$mu <- matrix(NA, nrow=init.b$K, ncol=dd$dimb)
        for (j in 1:dd$dimb) init.b$mu[,j] <- seq(bbmin[j]+afstand[j], bbmax[j]-afstand[j], length=init.b$K)
      }  
    }
    if (any(is.na(init.b$mu))) stop("NA in init.b$mu")          
    if (dd$dimb == 1){
      init.b$mu <- as.numeric(init.b$mu)
      if (length(init.b$mu) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$mu <- init.b$mu[1:init.b$K]          
      if (length(init.b$mu) != init.b$K) stop(paste("init.b$mu must be of length ", init.b$K, sep=""))
      names(init.b$mu) <- paste("mu", 1:init.b$K, sep="")    
    }else{
      if (!is.matrix(init.b$mu)) stop("init.b$mu must be a matrix")
      if (ncol(init.b$mu) != dd$dimb) stop(paste("init.b$mu must have ", dd$dimb, " columns", sep=""))
      if (nrow(init.b$mu) != init.b$K) stop(paste("init.b$mu must have ", init.b$K, " rows", sep=""))
      rownames(init.b$mu) <- paste("j", 1:init.b$K, sep="")
      colnames(init.b$mu) <- paste("m", 1:dd$dimb, sep="")        
    }

    ##### init.b:  Sigma and Li
    ##### ----------------------------------------------------
    if (dd$dimb == 1) bbVar <- (ifit$iSDranefVec / scale.b$scale)^2
    else              bbVar <- diag((ifit$iSDranefVec / scale.b$scale)^2)
    
    if (is.na(ib.Sigma)){            
      if (is.na(ib.Li)){       ### Sigma and Li are computed from the data
        if (dd$dimb == 1){
          init.b$Sigma <- rep(bbVar, init.b$K)
          names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")
          init.b$Li <- sqrt(1 / init.b$Sigma)
          names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
        }else{
          init.b$Sigma <- matrix(rep(t(bbVar), init.b$K), ncol=dd$dimb, byrow=TRUE)
          Sigmainv <- chol(bbVar)        
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init.b$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init.b$K)
          rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dd$dimb), ".", rep(1:dd$dimb, init.b$K), sep="")
          colnames(init.b$Sigma) <- paste("m", 1:dd$dimb, sep="")                
          names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=dd$LTb), rep(dd$naamLTb, init.b$K), sep="")
        }        
      }else{                 ### Li is checked and Sigma is computed from Li
        if (any(is.na(init.b$Li))) stop("NA in init.b$Li")                    
        if (dd$dimb == 1){
          if (length(init.b$Li) == 1) init.b$Li <- rep(init.b$Li, init.b$K)
          if (length(init.b$Li) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$Li <- init.b$Li[1:init.b$K]
          if (length(init.b$Sigma) != init.b$K) stop(paste("init.b$Sigma must be of length ", init.b$K, sep=""))
          init.b$Li <- as.numeric(init.b$Li)
          names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
          if (any(init.b$Li <= 0)) stop("init.b$Li must be positive")
          init.b$Sigma <- (1 / init.b$Li)^2
          names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")      
        }else{
          if (length(init.b$Li) == dd$LTb){
            tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
            tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init.b$Li
            tmpSigma <- tmpSigma %*% t(tmpSigma)
            err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
            if (class(err) == "try-error") stop("init.b$Li does not lead to a positive definite matrix")
            tmpSigma <- chol2inv(tmpSigma)
            init.b$Sigma <- matrix(rep(t(tmpSigma), init.b$K), ncol=dd$dimb, byrow=TRUE)
            init.b$Li <- rep(init.b$Li, init.b$K)
          }else{
            if (length(init.b$Li) == prior.b$Kmax*dd$LTb & prior.b$Kmax > init.b$K) init.b$Li <- init.b$Li[1:(init.b$K*dd$LTb)]
            if (length(init.b$Li) != init.b$K*dd$LTb) stop(paste("init.b$Li must be of length ", init.b$K*dd$LTb, sep=""))
            init.b$Sigma <- matrix(NA, ncol=dd$dimb, nrow=dd$dimb*init.b$K)
            for (j in 1:init.b$K){
              tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
              tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init.b$Li[((j-1)*dd$LTb+1):(j*dd$LTb)]
              tmpSigma <- tmpSigma %*% t(tmpSigma)
              err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
              if (class(err) == "try-error") stop(paste("the ", j,"-th block of init.b$Li does not lead to a positive definite matrix", sep=""))
              tmpSigma <- chol2inv(tmpSigma)
              init.b$Sigma[((j-1)*dd$dimb):(j*dd$dimb),] <- tmpSigma
            }  
          } 
          rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dd$dimb), ".", rep(1:dd$dimb, init.b$K), sep="")
          colnames(init.b$Sigma) <- paste("m", 1:dd$dimb, sep="")        
          names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=dd$LTb), rep(dd$naamLTb, init.b$K), sep="")
        }  
      }   
    }else{                   ### Sigma is checked and Li is computed from Sigma
      if (any(is.na(init.b$Sigma))) stop("NA in init.b$Sigma")              
      if (dd$dimb == 1){
        if (length(init.b$Sigma) == 1) init.b$Sigma <- rep(init.b$Sigma, init.b$K)
        if (length(init.b$Sigma) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$Sigma <- init.b$Sigma[1:init.b$K]      
        if (length(init.b$Sigma) != init.b$K) stop(paste("init.b$Sigma must be of length ", init.b$K, sep=""))
        init.b$Sigma <- as.numeric(init.b$Sigma)
        names(init.b$Sigma) <- paste("Sigma", 1:init.b$K, sep="")      
        if (any(init.b$Sigma <= 0)) stop("init.b$Sigma must be positive")
        init.b$Li <- sqrt(1 / init.b$Sigma)
        names(init.b$Li) <- paste("Li", 1:init.b$K, sep="")      
      }else{
        if (!is.matrix(init.b$Sigma)) stop("init.b$Sigma must be a matrix")
        if (ncol(init.b$Sigma) != dd$dimb) stop(paste("init.b$Sigma must have ", dd$dimb, " columns", sep=""))
        if (nrow(init.b$Sigma) == dd$dimb){
          if (any(init.b$Sigma[lower.tri(init.b$Sigma)] != t(init.b$Sigma)[lower.tri(init.b$Sigma)])) stop("init.b$Sigma must be a symmetric matrix")
          err <- try(Sigmainv <- chol(init.b$Sigma), silent=TRUE)
          if (class(err) == "try-error") stop("Cholesky decomposition of init.b$Sigma failed")
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init.b$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init.b$K)
        }else{
          if (nrow(init.b$Sigma) == prior.b$Kmax & prior.b$Kmax > init.b$K) init.b$Sigma <- init.b$Sigma[1:(init.b$K*dd$dimb),]
          if (nrow(init.b$Sigma) != init.b$K*dd$dimb) stop(paste("init.b$Sigma must have ", init.b$K, " times ", dd$dimb, " rows", sep=""))
          init.b$Li <- numeric(0)
          for (j in 1:init.b$K){
            Sigmainv <- init.b$Sigma[((j-1)*dd$dimb+1):(j*dd$dimb),]
            if (any(Sigmainv[lower.tri(Sigmainv)] != t(Sigmainv)[lower.tri(Sigmainv)])) stop(paste(j, "-th block of init.b$Sigma is not symmetric", sep=""))
            err <- try(Sigmainv <- chol(Sigmainv), silent=TRUE)
            if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of init.b$Sigma failed", sep=""))
            Sigmainv <- chol2inv(Sigmainv)
            Litmp <- t(chol(Sigmainv))
            init.b$Li <- c(init.b$Li, Litmp[lower.tri(Litmp, diag=TRUE)])
          }
        }
        rownames(init.b$Sigma) <- paste("j", rep(1:init.b$K, each=dd$dimb), ".", rep(1:dd$dimb, init.b$K), sep="")
        colnames(init.b$Sigma) <- paste("m", 1:dd$dimb, sep="")              
        names(init.b$Li) <- paste("Li", rep(1:init.b$K, each=dd$LTb), rep(dd$naamLTb, init.b$K), sep="")
      }  
    }    
    
    ##### init.b:  gammaInv
    ##### ----------------------------------------------------  
    if (is.na(ib.gammaInv)){
      if (dd$dimb == 1) init.b$gammaInv <- prior.b$zeta * bbVar
      else              init.b$gammaInv <- prior.b$zeta * diag(bbVar)
    }
    init.b$gammaInv <- as.numeric(init.b$gammaInv)
    if (length(init.b$gammaInv) == 1) init.b$gammaInv <- rep(init.b$gammaInv, dd$dimb)
    if (length(init.b$gammaInv) != dd$dimb) stop(paste("init.b$gammaInv must be of length ", dd$dimb, sep=""))
    if (any(is.na(init.b$gammaInv))) stop("NA in init.b$gammaInv")
    names(init.b$gammaInv) <- paste("gammaInv", 1:dd$dimb, sep="")
    
    ##### init.b:  r
    ##### ----------------------------------------------------  
    if (is.na(ib.r)){
      if (dd$dimb == 1){
        initz <- (init.b$b - scale.b$shift)/scale.b$scale
        MEANS <- matrix(rep(init.b$mu, ifit$I), ncol=init.b$K, byrow=TRUE)
        SDS   <- matrix(rep(sqrt(init.b$Sigma), ifit$I), ncol=init.b$K, byrow=TRUE)
        YY    <- matrix(rep(initz, init.b$K), ncol=init.b$K)
        WW    <- matrix(rep(init.b$w, ifit$I), ncol=init.b$K, byrow=TRUE)
        PROB  <- WW * dnorm(YY, mean=MEANS, sd=SDS)
      }else{
        initz <- (init.b$b - matrix(rep(scale.b$shift, ifit$I), ncol=dd$dimb, byrow=TRUE))/matrix(rep(scale.b$scale, ifit$I), ncol=dd$dimb, byrow=TRUE)
        PROB <- matrix(0, nrow=ifit$I, ncol=init.b$K)
        for (j in 1:init.b$K){
          MEANS <- init.b$mu[((j-1)*dd$dimb+1):(j*dd$dimb)]
          SIGMA <- init.b$Sigma[((j-1)*dd$dimb+1):(j*dd$dimb),]
          PROB[,j] <- init.b$w[j] * dMVN(initz, mean=MEANS, Sigma=SIGMA)        
        }        
      }
      sumPROB <- apply(PROB, 1, sum)
      sumPROB[sumPROB <= 0] <- 1
      PROB    <- PROB / matrix(rep(sumPROB, each=init.b$K), ncol=init.b$K, byrow=TRUE)
      init.b$r <- apply(PROB, 1, which.max)          
    }
    init.b$r <- as.numeric(init.b$r)
    if (length(init.b$r) != ifit$I) stop(paste("init.b$r must be of length ", ifit$I, sep=""))
    if (any(is.na(init.b$r))) stop("NA in init.b$r")
    if (any(init.b$r < 1) | any(init.b$r > init.b$K)) stop(paste("init.b$r out of the range (must lie between ", 1, " and ", init.b$K, ")", sep=""))
    names(init.b$r) <- paste("r", 1:ifit$I, sep="")    
  }else{
    init.b <- list(b=0, K=0, w=0, mu=0, Sigma=0, Li=0, gammaInv=0, r=0)
  }  


########## ========== Parameters from inits (to allow for vectorized calculation in a future) ========== ##########
########## ============================================================================================= ##########
  Csigma_eps <- init.eps$sigma
  CgammaInv_eps <- init.eps$gammaInv
  if (dd$dimb){  
    CK_b <- init.b$K
    Cw_b <- c(init.b$w, rep(0, prior.b$Kmax - init.b$K))
    if (dd$dimb == 1){
      Cmu_b<- c(init.b$mu, rep(0, prior.b$Kmax - init.b$K))
      CLi_b<- c(init.b$Li, rep(0, prior.b$Kmax - init.b$K))
    
    }else{
      Cmu_b<- c(t(init.b$mu), rep(0, dd$dimb*(prior.b$Kmax - init.b$K)))
      CLi_b<- c(init.b$Li, rep(0, dd$LTb*(prior.b$Kmax - init.b$K)))
    }
    CgammaInv_b <- init.b$gammaInv
    Cr_b <- init.b$r - 1    
    Cbb  <- as.numeric(t(init.b$b))
  }else{
    CK_b        <- 0
    Cw_b        <- 0
    Cmu_b       <- 0
    CLi_b       <- 0
    CgammaInv_b <- 0
    Cr_b        <- 0
    Cbb         <- 0
  }  
  Cbeta <- init.beta


########## ========== nMCMC ========== ##########
########## =========================== ##########
  if (length(nMCMC) != 4) stop("nMCMC must be of length 4")
  if (is.null(names(nMCMC))) names(nMCMC) <- c("burn", "keep", "thin", "info")
  names.nMCMC <- names(nMCMC)
  if (!match("burn", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("burn"), "] must be specified", sep=""))
  else                                        n.burn <- nMCMC["burn"]
  if (!match("keep", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("keep"), "] must be specified", sep=""))
  else                                        n.keep <- nMCMC["keep"]
  if (!match("thin", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("thin"), "] must be specified", sep=""))
  else                                        n.thin <- nMCMC["thin"]
  if (!match("info", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("info"), "] must be specified", sep=""))
  else                                        n.info <- nMCMC["info"]
  nMCMC <- c(n.burn, n.keep, n.thin, n.info)
  names(nMCMC) <- c("burn", "keep", "thin", "info")  
  if (nMCMC["burn"] < 0) stop(paste("nMCMC[", dQuote("burn"), "] must be non-negative", sep=""))
  if (nMCMC["keep"] <= 0) stop(paste("nMCMC[", dQuote("keep"), "] must be positive", sep=""))  
  if (nMCMC["thin"] <= 0) stop(paste("nMCMC[", dQuote("thin"), "] must be positive", sep=""))
  if (nMCMC["info"] <= 0 | nMCMC["info"] > max(nMCMC["burn"], nMCMC["keep"])) nMCMC["info"] <- max(nMCMC["burn"], nMCMC["keep"])


########## ========== store ========== ##########
########## =========================== ##########
  if (length(store) != 1) stop("store must be of length 1")
  if (is.null(names(store))) names(store) <- c("b")
  names.store <- names(store)
  if (!match("b", names.store, nomatch=0)) stop(paste("store[", dQuote("b"), "] must be specified", sep=""))
  else                                     store.b <- store["b"]
  store <- c(store.b)
  names(store) <- c("b")
  if (!dd$dimb) store["b"] <- FALSE  
  

########## ========== Show values of some variables when debugging ========== ##########
########## ================================================================== ##########
  DEBUG <- FALSE
  if (DEBUG){
    cat("ifit$n:\n------------\n")
    print(ifit$n)
    
    cat("\nifit$Cn:\n------------\n")
    print(ifit$Cn)

    cat("\ndd$dimb: ", dd$dimb, "\n-----------\n")

    cat("\ndd$lbeta: ", dd$lbeta, "\n-----------\n")    

    cat("\ndd$p:\n-----------\n")
    print(dd$p)

    cat("\ndd$q:\n-----------\n")
    print(dd$q)    

    cat("\ndd$random.intercept:\n-----------\n")
    print(dd$random.intercept)    

    cat("\ndd$CrandomIntcpt:\n-----------\n")
    print(dd$CrandomIntcpt)    

    cat("\ndd$fixed.intercept:\n-----------\n")
    print(dd$fixed.intercept)    

    cat("\ndd$CfixedIntcpt:\n-----------\n")
    print(dd$CfixedIntcpt)    
    
    cat("\nifit$Cy_c:\n------------\n")
    print(ifit$Cy_c)

    cat("\nifit$Cy_d:\n------------\n")
    print(ifit$Cy_d)
    
    cat("\nifit$CX:\n------------\n")
    print(ifit$CX)

    cat("\nifit$CZ:\n------------\n")
    print(ifit$CZ)

    cat("\nifit$CXtX:\n-----------\n")
    print(ifit$CXtX)
    
    cat("\nifit$CZitZi:\n-----------\n")
    print(ifit$CZitZi)
    
    cat("\nifit$iintcpt:\n------------\n")
    print(ifit$iintcpt)

    cat("\nifit$is.intcpt:\n------------\n")
    print(ifit$is.intcpt)

    cat("\nifit$ifixef:\n------------\n")
    print(ifit$ifixef)

    cat("ifit$is.fixef:\n------------\n")
    print(ifit$is.fixef)

    cat("\nifit$ibeta:\n-----------\n")
    print(ifit$ibeta)
    
    cat("\nifit$isigma:\n------------\n")
    print(ifit$isigma)

    cat("\nifit$is.sigma:\n------------\n")
    print(ifit$is.sigma)

    cat("ifit$is.ranef:\n------------\n")
    print(ifit$is.ranef)
    
    cat("\nifit$iEranefVec:\n------------\n")
    print(ifit$iEranefVec)
    
    cat("\nifit$iEranef:\n------------\n")
    print(ifit$iEranef)

    cat("\nifit$iSEranefVec:\n------------\n")
    print(ifit$iSEranefVec)
    
    cat("\nifit$iSDranefVec:\n------------\n")
    print(ifit$iSDranefVec)
    
    cat("\nifit$iSDranef:\n------------\n")
    print(ifit$iSDranef)

    cat("\nifit$iSEranefVec:\n------------\n")
    print(ifit$iSEranefVec)

    cat("\nifit$ibMat:\n------------\n")
    print(ifit$ibMat)
    
    cat("\nscale.b:\n------------\n")
    print(scale.b)
    
    cat("scb$CshiftScale_b:\n------------\n")
    print(scb$CshiftScale_b)

    cat("pbeta$CpriorDouble_beta:\n-----------\n")
    print(pbeta$CpriorDouble_beta)
    
    cat("peps$CpriorDouble_eps:\n-----------\n")
    print(peps$CpriorDouble_eps)

    cat("pbb$CpriorInt_b:\n-----------\n")
    print(pbb$CpriorInt_b)    
    cat("pbb$CpriorDouble_b:\n-----------\n")
    print(pbb$CpriorDouble_b)

    cat("Csigma_eps:\n-----------\n")
    print(Csigma_eps)    
    cat("CgammaInv_eps:\n-----------\n")
    print(CgammaInv_eps)    

    cat("CK_b:\n-----------\n")
    print(CK_b)    
    cat("Cw_b:\n-----------\n")
    print(Cw_b)    
    cat("Cmu_b:\n-----------\n")
    print(Cmu_b)    
    cat("CLi_b:\n-----------\n")
    print(CLi_b)            
    cat("CgammaInv_b:\n-----------\n")
    print(CgammaInv_b)    
    cat("Cr_b:\n-----------\n")
    print(Cr_b)    
    cat("Cbb:\n-----------\n")
    print(Cbb)    
    
    cat("\n")
  }  


  SHOW.CLUST <- FALSE
  if (SHOW.CLUST){
    cumq_ri <- cumsum(dd$q_ri)
    ishow <- 98
    IDshow <- unique(dd$id)[ishow]
    
    iY <- iX <- iZ <- iZstar <- efixed <- erand <- ezs <- list()
    inn <- numeric(dd$R)
    for (s in 1:dd$R){
      iRow <- (ifit$ID[[s]] == IDshow)
      if (s == 1) iRand <- if (dd$q_ri[s]) 1:dd$q_ri[s] else NULL
      else        iRand <- if (dd$q_ri[s]) (cumq_ri[s-1]+1):cumq_ri[s] else NULL
      inn[s] <- ifit$n[[s]][ishow]
      
      iY[[s]] <- ifit$Y[[s]][iRow]
      
      if (!is.character(ifit$x[[s]])){
        if (ncol(ifit$x[[s]]) == 1){
          efixed[[s]] <- ifit$x[[s]][iRow,] * as.numeric(ifit$ifixef[[s]][, "Est"])
          iifit$x[[s]] <- matrix(ifit$x[[s]][iRow,], ncol=1)
        }else{  
          efixed[[s]] <- as.numeric(ifit$x[[s]][iRow,] %*% as.numeric(ifit$ifixef[[s]][, "Est"]))
          iX[[s]] <- ifit$x[[s]][iRow,]
        }  
      }else{
        efixed[[s]] <- rep(0, sum(iRow))
        iX[[s]] <- "empty"
      }  
      if (!is.character(ifit$z[[s]])){
        if (ncol(ifit$z[[s]]) == 1){
          erand[[s]] <- ifit$z[[s]][iRow,] * as.numeric(ifit$ib[[s]][as.character(IDshow),])
          ezs[[s]]   <- ifit$z[[s]][iRow,] * scale.b$scale[iRand]
          iZ[[s]] <- matrix(ifit$z[[s]][iRow,], ncol=1)          
        }else{  
          erand[[s]] <- as.numeric(ifit$z[[s]][iRow,] %*% as.numeric(ifit$ib[[s]][as.character(IDshow),]))
          ezs[[s]]   <- as.numeric(ifit$z[[s]][iRow,] %*% scale.b$scale[iRand])
          iZ[[s]] <- matrix(ifit$z[[s]][iRow,], ncol=ncol(ifit$z[[s]]))          
        }
        if (dd$R == 1){
          iZstar[[s]] <- iZ[[s]]
        }else{
          if (s == 1){
            if (sum(dd$q_ri[2:dd$R])) iZstar[[s]] <- cbind(iZ[[s]], matrix(0, ncol=sum(dd$q_ri[2:dd$R]), nrow=inn[s])) else iZstar[[s]] <- iZ[[s]]
          }else{
            if (s == dd$R){
              if (sum(dd$q_ri[1:(dd$R-1)])) iZstar[[s]] <- cbind(matrix(0, ncol=sum(dd$q_ri[1:(dd$R-1)]), nrow=inn[s]), iZ[[s]]) else iZstar[[s]] <- iZ[[s]]
            }else{
              if (sum(dd$q_ri[1:(s-1)])) iZstar[[s]] <- cbind(matrix(0, ncol=sum(dd$q_ri[1:(s-1)]), nrow=inn[s]), iZ[[s]]) else iZstar[[s]] <- iZ[[s]]
              if (sum(dd$q_ri[(s+1):dd$R])) iZstar[[s]] <- cbind(iZstar[[s]], matrix(0, ncol=sum(dd$q_ri[(s+1):dd$R]), nrow=inn[s])) else iZstar[[s]] <- iZstar[[s]]
            }  
          }  
        }  
      }else{
        erand[[s]] <- ezs[[s]] <- rep(0, sum(iRow))
        iZ[[s]] <- "empty"
        iZstar[[s]] <- matrix(0, ncol=sum(dd$q_ri), nrow=inn[s])
      }  
    }  

    iZZstar <- iZstar[[1]]
    if (dd$R > 1) for (s in 2:dd$R) iZZstar <- rbind(iZZstar, iZstar[[s]])
    
    MU <- init.b$mu[init.b$r[ishow],]
    SIGMA <- init.b$Sigma[((init.b$r[ishow]-1)*dd$dimb+1):(init.b$r[ishow]*dd$dimb),]
    Q <- chol2inv(chol(SIGMA))

    if (dd$dimb > 1) Smat <- diag(scale.b$scale) else Smat <- matrix(scale.b$scale, nrow=1, ncol=1)

    cat("Cluster ", ishow, ": \n", sep="")
    for (s in 1:dd$R){
      cat("s=", s, ", SZ'ZS:\n", sep="")
      if (dd$q_ri[s]){
        if (s == 1) indr <- 1:dd$q_ri[s] else indr <- (cumq_ri[s-1]+1):cumq_ri[s]
        Ss <- matrix(Smat[indr, indr], nrow=dd$q_ri[s], ncol=dd$q_ri[s])
        print(Ss %*% t(iZ[[s]]) %*% iZ[[s]] %*% Ss)
      }else{
        cat("empty\n")
      }  
    }      
  }  

  
########## ========== MCMC simulation                              ========== ##########
########## ================================================================== ##########
  cat(paste("MCMC sampling started on ", date(), ".\n", sep=""))  
  MCMC <- .C("GLMM_MCMC",
             Y_c              = as.double(ifit$Cy_c),
             R_c              = as.integer(dd$Rc),
             Y_d              = as.integer(ifit$Cy_d),
             R_d              = as.integer(dd$Rd),
             dist             = as.integer(dd$ndist),
             I                = as.integer(ifit$I),
             n                = as.integer(ifit$Cn),
             X                = as.double(ifit$CX),
             XtX              = as.double(ifit$CXtX),
             p                = as.integer(dd$p),
             fixedIntcpt      = as.integer(dd$CfixedIntcpt),
             Z                = as.double(ifit$CZ),
             ZitZi            = as.double(ifit$CZitZi),
             q                = as.integer(dd$q),
             randIntcpt       = as.integer(dd$CrandomIntcpt),
             shiftScale_b     = as.double(scb$CshiftScale_b),
             nMCMC            = as.integer(nMCMC),
             keepChain        = as.integer(store),
             priorDouble_eps  = as.double(peps$CpriorDouble_eps),
             priorInt_b       = as.integer(pbb$CpriorInt_b),
             priorDouble_b    = as.double(pbb$CpriorDouble_b),
             priorDouble_beta = as.double(pbeta$CpriorDouble_beta),
             sigma_eps        = as.double(Csigma_eps),
             gammaInv_eps     = as.double(CgammaInv_eps),
             K_b              = as.integer(CK_b),
             w_b              = as.double(Cw_b),
             mu_b             = as.double(Cmu_b),
             Q_b              = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             Sigma_b          = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             Li_b             = as.double(CLi_b),
             gammaInv_b       = as.double(CgammaInv_b),
             r_b              = as.integer(Cr_b),
             beta             = as.double(Cbeta),
             b                = as.double(Cbb),
             chsigma_eps      = double(ifelse(dd$Rc, dd$Rc * nMCMC["keep"], 1)),
             chgammaInv_eps   = double(ifelse(dd$Rc, dd$Rc * nMCMC["keep"], 1)),
             chK_b            = integer(ifelse(dd$dimb, nMCMC["keep"], 1)),
             chw_b            = double(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chmu_b           = double(ifelse(dd$dimb, dd$dimb * prior.b$Kmax * nMCMC["keep"], 1)),
             chQ_b            = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chSigma_b        = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chLi_b           = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chgammaInv_b     = double(ifelse(dd$dimb, dd$dimb * nMCMC["keep"], 1)),
             chorder_b        = integer(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chrank_b         = integer(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chMeanData_b     = double(ifelse(dd$dimb, dd$dimb * nMCMC["keep"], 1)),
             chCorrData_b     = double(ifelse(dd$dimb, dd$LTb * nMCMC["keep"], 1)),
             chbeta           = double(ifelse(dd$lbeta, dd$lbeta * nMCMC["keep"], 1)),
             chb              = double(ifelse(dd$dimb, ifelse(store["b"], ifit$I * dd$dimb * nMCMC["keep"], ifit$I * dd$dimb), 1)),
             pm_eta_fixed     = double(ifit$sumCn),
             pm_eta_random    = double(ifit$sumCn),             
             pm_b             = double(ifelse(dd$dimb, dd$dimb * ifit$I, 1)),
             pm_w_b           = double(ifelse(dd$dimb, prior.b$Kmax, 1)),
             pm_mu_b          = double(ifelse(dd$dimb, dd$dimb * prior.b$Kmax, 1)),
             pm_Q_b           = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_Sigma_b       = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_Li_b          = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_indLogL       = double(ifit$I),
             pm_indLogpb      = double(ifelse(dd$dimb, ifit$I, 1)),
             iter             = as.integer(0),
             err              = as.integer(0),
             PACKAGE=thispackage)            
  cat(paste("MCMC sampling finished on ", date(), ".\n", sep=""))
  if (MCMC$err) stop("Something went wrong.")


  ########## ========== State of MCMC ========== ##########
  ########## =================================== ##########
  if (dd$dimb){
    state.w_b <- as.numeric(MCMC$w_b[1:MCMC$K_b])
    names(state.w_b) <- paste("w", 1:MCMC$K_b, sep="")
    state.r_b <- as.numeric(MCMC$r_b + 1)
    names(state.r_b) <- paste("r", 1:ifit$I, sep="")
    state.gammaInv_b <- as.numeric(MCMC$gammaInv_b)
    names(state.gammaInv_b) <- paste("gammaInv", 1:dd$dimb, sep="")
    if (dd$dimb == 1){
      state.mu_b <- as.numeric(MCMC$mu_b[1:MCMC$K_b])
      names(state.mu_b) <- paste("mu", 1:MCMC$K_b, sep="")
      state.Li_b <- as.numeric(MCMC$Li_b[1:MCMC$K_b])
      names(state.Li_b) <- paste("Li", 1:MCMC$K_b, sep="")
      state.Sigma_b <- (1 / state.Li_b)^2
      names(state.Sigma_b) <- paste("Sigma", 1:MCMC$K_b, sep="")
      state.Q_b <- as.numeric(MCMC$Q_b[1:MCMC$K_b])
      names(state.Q_b) <- paste("Q", 1:MCMC$K_b, sep="")

      state.b <- as.numeric(MCMC$b)
      names(state.b) <- 1:ifit$I
    }else{
      state.mu_b <- matrix(MCMC$mu_b[1:(dd$dimb*MCMC$K_b)], ncol=dd$dimb, byrow=TRUE)
      rownames(state.mu_b) <- paste("j", 1:MCMC$K_b, sep="")
      colnames(state.mu_b) <- paste("m", 1:dd$dimb, sep="")
      state.Li_b <- as.numeric(MCMC$Li_b[1:(dd$LTb*MCMC$K_b)])
      names(state.Li_b) <- paste("Li", rep(1:MCMC$K_b, each=dd$LTb), rep(dd$naamLTb, MCMC$K_b), sep="")
      state.Sigma_b<- matrix(NA, ncol=dd$dimb, nrow=dd$dimb*MCMC$K_b)
      rownames(state.Sigma_b) <- paste("j", rep(1:MCMC$K_b, each=dd$dimb), ".", rep(1:dd$dimb, MCMC$K_b), sep="")
      colnames(state.Sigma_b) <- paste("m", 1:dd$dimb, sep="")            
      for (j in 1:MCMC$K_b){
        tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state.Li_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpSigma <- tmpSigma %*% t(tmpSigma)
        tmpSigma <- chol2inv(chol(tmpSigma))
        state.Sigma_b[((j-1)*dd$dimb+1):(j*dd$dimb),] <- tmpSigma
      }
      state.Q_b<- as.numeric(MCMC$Q_b[1:(dd$LTb*MCMC$K_b)])
      names(state.Q_b) <- paste("Q", rep(1:MCMC$K_b, each=dd$LTb), rep(dd$naamLTb, MCMC$K_b), sep="")

      state.b <- matrix(MCMC$b, ncol=dd$dimb, nrow=ifit$I, byrow=TRUE)
      colnames(state.b) <- paste("b", 1:dd$dimb, sep="")
      rownames(state.b) <- 1:ifit$I
      
    }
    nCompTotal_b<- sum(MCMC$chK_b)
    freqK_b <- table(MCMC$chK_b)
    propK_b <- prop.table(freqK_b)
  }else{
    state.w_b <- state.r_b <- state.gamma_b <- state.mu_b <- state.Li_b <- state.Sigma_b <- state.Q_b <- state.b <- 0
  }  

  if (dd$lbeta){
    state.beta <- as.numeric(MCMC$beta)
    names(state.beta) <- paste("beta", 1:dd$lbeta, sep="")
  }else{
    state.beta <- 0
  }  

  if (dd$Rc){
    state.sigma_eps <- as.numeric(MCMC$sigma_eps)
    names(state.sigma_eps) <- paste("sigma", 1:dd$Rc, sep="")
    state.gammaInv_eps <- as.numeric(MCMC$gammaInv_eps)
    names(state.gammaInv_eps) <- paste("gammaInv", 1:dd$Rc, sep="")
  }else{
    state.sigma_eps <- state.gammaInv_eps <- 0
  }  

  
  ########## ========== Create a list to be returned ========== ##########
  ########## ================================================== ##########
  RET <- list(iter             = MCMC$iter,
              nMCMC            = nMCMC,
              dist             = dd$dist,
              R                = c(Rc=dd$Rc, Rd=dd$Rd),
              p                = dd$p,
              q                = dd$q,
              fixed.intercept  = dd$fixed.intercept,
              random.intercept = dd$random.intercept,
              lbeta            = dd$lbeta,
              dimb             = dd$dimb,
              prior.beta       = prior.beta,
              prior.b          = prior.b,
              prior.eps        = prior.eps,
              init.beta        = init.beta,
              init.b           = init.b,
              init.eps         = init.eps,
              state.beta       = state.beta,
              state.b          = list(b        = state.b,
                                      K        = as.numeric(MCMC$K_b),
                                      w        = state.w_b,
                                      mu       = state.mu_b,
                                      Sigma    = state.Sigma_b,
                                      Li       = state.Li_b,
                                      Q        = state.Q_b,
                                      gammaInv = state.gammaInv_b,
                                      r        = state.r_b),
              state.eps        = list(sigma    = state.sigma_eps,
                                      gammaInv = state.gammaInv_eps),
              scale.b          = scale.b,
              freqK_b          = freqK_b,
              propK_b          = propK_b)

  
  ########## ========== Posterior means of quantities computed in C++ ========== ##########
  ########## =================================================================== ##########
  RET$poster.mean.eta <- data.frame(fixed  = as.numeric(MCMC$pm_eta_fixed),
                                    random = as.numeric(MCMC$pm_eta_random))

  if (dd$dimb){
    MCMC$pm_b <- matrix(MCMC$pm_b, ncol=dd$dimb, byrow=TRUE)
    RET$poster.mean.cluster <- as.data.frame(MCMC$pm_b)
    colnames(RET$poster.mean.cluster) <- paste("b", 1:dd$dimb, sep="")    
    RET$poster.mean.cluster$LogL  <- as.numeric(MCMC$pm_indLogL)
    RET$poster.mean.cluster$Logpb <- as.numeric(MCMC$pm_indLogpb)
    
    if (prior.b$priorK == "fixed"){
                                                      ##### I am not sure whether the posterior means (especially of variance components) are useful!
                                                      ##### In any case, they should be used with care
                                                      ##### -----------------------------------------------------------------------------------------
      RET$poster.mean.w_b <- as.numeric(MCMC$pm_w_b)
      names(RET$poster.mean.w_b) <- paste("w", 1:prior.b$Kmax, sep="")

      RET$poster.mean.mu_b <- matrix(MCMC$pm_mu_b, nrow=prior.b$Kmax, ncol=dd$dimb, byrow=TRUE)
      rownames(RET$poster.mean.mu_b) <- paste("j", 1:prior.b$Kmax, sep="")
      colnames(RET$poster.mean.mu_b) <- paste("m", 1:dd$dimb, sep="")

      RET$poster.mean.Q_b <- RET$poster.mean.Sigma_b <- RET$poster.mean.Li_b <- list()
      for (j in 1:prior.b$Kmax){
        tmpQ <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpQ[lower.tri(tmpQ, diag=TRUE)] <- MCMC$pm_Q_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpQ[upper.tri(tmpQ, diag=FALSE)] <- t(tmpQ)[upper.tri(t(tmpQ), diag=FALSE)]
        RET$poster.mean.Q_b[[j]] <- tmpQ
      
        tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- MCMC$pm_Sigma_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpSigma[upper.tri(tmpSigma, diag=FALSE)] <- t(tmpSigma)[upper.tri(t(tmpSigma), diag=FALSE)]
        RET$poster.mean.Sigma_b[[j]] <- tmpSigma
      
        tmpLi <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpLi[lower.tri(tmpLi, diag=TRUE)] <- MCMC$pm_Li_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpLi[upper.tri(tmpLi, diag=FALSE)] <- t(tmpLi)[upper.tri(t(tmpLi), diag=FALSE)]
        RET$poster.mean.Li_b[[j]] <- tmpLi      
      }
      names(RET$poster.mean.Q_b) <- names(RET$poster.mean.Sigma_b) <- names(RET$poster.mean.Li_b) <- paste("j", 1:prior.b$Kmax, sep="")    
    }      
  }else{
    RET$poster.mean.cluster <- data.frame(LogL  = as.numeric(MCMC$pm_indLogL))
  }  

  ########## ========== Additional posterior summaries                ========== ##########
  ########## =================================================================== ##########
  qProbs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  nSumm <- c("Mean", "Std.Dev.", "Min.", "2.5%", "1st Qu.", "Median", "3rd Qu.", "97.5%", "Max.")

  if (dd$lbeta){
    MCMC$chbeta <- matrix(MCMC$chbeta, ncol=dd$lbeta, byrow=TRUE)
    colnames(MCMC$chbeta) <- paste("beta", 1:dd$lbeta, sep="")
    
    if (dd$lbeta == 1){
      mean.beta  <- mean(MCMC$chbeta, na.rm=TRUE)
      quant.beta <- quantile(MCMC$chbeta, prob=qProbs, na.rm=TRUE)
      sd.beta    <- sd(MCMC$chbeta, na.rm=TRUE)
      RET$summ.beta     <- c(mean.beta, sd.beta, quant.beta)
      names(RET$summ.beta) <- nSumm      
    }else{
      mean.beta  <- apply(MCMC$chbeta, 2, mean, na.rm=TRUE)
      quant.beta <- apply(MCMC$chbeta, 2, quantile, prob=qProbs, na.rm=TRUE)
      sd.beta    <- apply(MCMC$chbeta, 2, sd, na.rm=TRUE)
      RET$summ.beta     <- rbind(mean.beta, sd.beta, quant.beta)
      rownames(RET$summ.beta) <- nSumm            
    }
  }  

  if (dd$dimb){
    MCMC$chMeanData_b <- matrix(MCMC$chMeanData_b, ncol=dd$dimb, byrow=TRUE)
    MCMC$chCorrData_b <- matrix(MCMC$chCorrData_b, ncol=dd$LTb, byrow=TRUE)
    colnames(MCMC$chMeanData_b) <- paste("b.Mean.", 1:dd$dimb, sep="")
    colnames(MCMC$chCorrData_b) <- paste("b.Corr", dd$naamLTb, sep="")
    colnames(MCMC$chCorrData_b)[((0:(dd$dimb-1))*(2*dd$dimb - (0:(dd$dimb-1)) + 1))/2 + 1] <- paste("b.SD.", 1:dd$dimb, sep="")  
    
    if (dd$dimb == 1){
      meanb.Mean  <- mean(MCMC$chMeanData_b, na.rm=TRUE)
      quantb.Mean <- quantile(MCMC$chMeanData_b, prob=qProbs, na.rm=TRUE)    
      sdb.Mean    <- sd(MCMC$chMeanData_b, na.rm=TRUE)
      RET$summ.b.Mean <- c(meanb.Mean, sdb.Mean, quantb.Mean)
      names(RET$summ.b.Mean) <- nSumm

      meanb.SDCorr  <- mean(MCMC$chCorrData_b, na.rm=TRUE)
      quantb.SDCorr <- quantile(MCMC$chCorrData_b, prob=qProbs, na.rm=TRUE)    
      sdb.SDCorr    <- sd(MCMC$chCorrData_b, na.rm=TRUE)
      RET$summ.b.SDCorr <- c(meanb.SDCorr, sdb.SDCorr, quantb.SDCorr)
      names(RET$summ.b.SDCorr) <- nSumm    
    }else{
      meanb.Mean  <- apply(MCMC$chMeanData_b, 2, mean, na.rm=TRUE)
      quantb.Mean <- apply(MCMC$chMeanData_b, 2, quantile, prob=qProbs, na.rm=TRUE)    
      sdb.Mean    <- apply(MCMC$chMeanData_b, 2, sd, na.rm=TRUE)
      RET$summ.b.Mean <- rbind(meanb.Mean, sdb.Mean, quantb.Mean)
      rownames(RET$summ.b.Mean) <- nSumm

      meanb.SDCorr  <- apply(MCMC$chCorrData_b, 2, mean, na.rm=TRUE)
      quantb.SDCorr <- apply(MCMC$chCorrData_b, 2, quantile, prob=qProbs, na.rm=TRUE)    
      sdb.SDCorr    <- apply(MCMC$chCorrData_b, 2, sd, na.rm=TRUE)
      RET$summ.b.SDCorr <- rbind(meanb.SDCorr, sdb.SDCorr, quantb.SDCorr)
      rownames(RET$summ.b.SDCorr) <- nSumm        
    }  
  }  

  if (dd$Rc){
    MCMC$chsigma_eps <- matrix(MCMC$chsigma_eps, ncol=dd$Rc, byrow=TRUE)
    colnames(MCMC$chsigma_eps) <- paste("sigma", 1:dd$Rc, sep="")

    if (dd$Rc == 1){
      mean.sigma_eps  <- mean(MCMC$chsigma_eps, na.rm=TRUE)
      quant.sigma_eps <- quantile(MCMC$chsigma_eps, prob=qProbs, na.rm=TRUE)
      sd.sigma_eps    <- sd(MCMC$chsigma_eps, na.rm=TRUE)
      RET$summ.sigma_eps     <- c(mean.sigma_eps, sd.sigma_eps, quant.sigma_eps)
      names(RET$summ.sigma_eps) <- nSumm      
    }else{
      mean.sigma_eps  <- apply(MCMC$chsigma_eps, 2, mean, na.rm=TRUE)
      quant.sigma_eps <- apply(MCMC$chsigma_eps, 2, quantile, prob=qProbs, na.rm=TRUE)
      sd.sigma_eps    <- apply(MCMC$chsigma_eps, 2, sd, na.rm=TRUE)
      RET$summ.sigma_eps     <- rbind(mean.sigma_eps, sd.sigma_eps, quant.sigma_eps)
      rownames(RET$summ.sigma_eps) <- nSumm            
    }    
  }  
  
  
  ########## ========== Chains for model parameters ========== ##########
  ########## ================================================= ##########
  if (keep.chains){
    if (dd$dimb){
      ##### Chains for parameters of mixture distribution of b
      ##### -----------------------------------------------------
      RET$K_b <- as.numeric(MCMC$chK_b)
      MCMC$K_b <- NULL

      RET$w_b <- as.numeric(MCMC$chw_b[1:nCompTotal_b])
      MCMC$chw_b <- NULL

      RET$mu_b <- as.numeric(MCMC$chmu_b[1:(dd$dimb*nCompTotal_b)])
      MCMC$chmu_b <- NULL

      RET$Li_b <- as.numeric(MCMC$chLi_b[1:(dd$LTb*nCompTotal_b)])
      MCMC$chLi_b <- NULL    

      RET$Q_b <- as.numeric(MCMC$chQ_b[1:(dd$LTb*nCompTotal_b)])
      MCMC$chQ_b <- NULL

      RET$Sigma_b <- as.numeric(MCMC$chSigma_b[1:(dd$LTb*nCompTotal_b)])
      MCMC$chSigma_b <- NULL

      RET$gammaInv_b <- matrix(MCMC$chgammaInv_b, ncol=dd$dimb, byrow=TRUE)
      colnames(RET$gammaInv_b) <- paste("gammaInv", 1:dd$dimb, sep="")  
      MCMC$chgammaInv_b <- NULL
  
      RET$order_b <- as.numeric(MCMC$chorder_b[1:nCompTotal_b] + 1)
      MCMC$chorder_b <- NULL

      RET$rank_b <- as.numeric(MCMC$chrank_b[1:nCompTotal_b] + 1)
      MCMC$chrank_b <- NULL

      if (prior.b$priorK == "fixed"){
        RET$w_b <- matrix(RET$w_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$w_b) <- paste("w", 1:prior.b$Kmax, sep="")

        RET$mu_b <- matrix(RET$mu_b, ncol=dd$dimb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$mu_b) <- paste("mu.", rep(1:prior.b$Kmax, each=dd$dimb), ".", rep(1:dd$dimb, prior.b$Kmax), sep="")
    
        RET$Li_b <- matrix(RET$Li_b, ncol=dd$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Li_b) <- paste("Li", rep(1:prior.b$Kmax, each=dd$LTb), rep(dd$naamLTb, prior.b$Kmax), sep="")

        RET$Q_b <- matrix(RET$Q_b, ncol=dd$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Q_b) <- paste("Q", rep(1:prior.b$Kmax, each=dd$LTb), rep(dd$naamLTb, prior.b$Kmax), sep="")

        RET$Sigma_b <- matrix(RET$Sigma_b, ncol=dd$LTb*prior.b$Kmax, byrow=TRUE)
        colnames(RET$Sigma_b) <- paste("Sigma", rep(1:prior.b$Kmax, each=dd$LTb), rep(dd$naamLTb, prior.b$Kmax), sep="")

        RET$order_b <- matrix(RET$order_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$order_b) <- paste("order", 1:prior.b$Kmax, sep="")

        RET$rank_b <- matrix(RET$rank_b, ncol=prior.b$Kmax, byrow=TRUE)
        colnames(RET$rank_b) <- paste("rank", 1:prior.b$Kmax, sep="")        
      }
            
      ##### Chains for characteristics of the mixture distribution of b
      ##### --------------------------------------------------------------
      RET$mixture_b <- as.data.frame(cbind(MCMC$chMeanData_b, MCMC$chCorrData_b))
      MCMC$chMeanData_b <- NULL
      MCMC$chCorrData_b <- NULL
      
      ##### Chains for random effects b
      ##### ------------------------------
      if (store.b){
        RET$b <- matrix(MCMC$b, ncol=dd$dimb*ifit$I, byrow=TRUE)
        MCMC$chb <- NULL
        colnames(RET$b) <- paste("b.", rep(1:ifit$I, each=dd$dimb), ".", rep(1:dd$dimb, ifit$I), sep="")
      }        
    }

    if (dd$lbeta){
      ##### Chains for regression coefficients beta
      ##### ------------------------------------------
      RET$beta <- MCMC$chbeta
      MCMC$chbeta <- NULL      
    }

    if (dd$Rc){
      ##### Chains for parameters of distribution of residuals
      ##### -----------------------------------------------------
      RET$sigma_eps <- MCMC$chsigma_eps
      MCMC$chsigma_eps <- NULL

      RET$gammaInv_eps <- matrix(MCMC$chgammaInv_eps, ncol=dd$Rc, byrow=TRUE)
      colnames(RET$gammaInv_eps) <- paste("gammaInv", 1:dd$Rc, sep="")
      MCMC$chgammaInv_eps <- NULL
    }        
  }    
  
  class(RET) <- "GLMM_MCMC"
  return(RET)    
}

