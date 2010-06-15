##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of continuous response and random effects
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    06/07/2009
##              03/08/2009:  version for continuous responses working
##              10/11/2009:  version for combined discrete and continuous responses working
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
                      tuneMCMC=list(beta=1, b=1),
                      store=c(b=FALSE), keep.chains=TRUE)
{
  require("lme4")
  thispackage <- "mixAK"

  DEBUG <- FALSE
  
  parallel <- FALSE     ### at this moment, parallel computation not fully implemented

########## ========== Data ========== ##########
########## ========================== ##########
  dd <- GLMM_MCMCdata(y=y, dist=dist, id=id, x=x, z=z, random.intercept=random.intercept)
  if (dd$R > 1) NAAM.RESP <- names(x)
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
          if (nrow(init.b$Sigma) == prior.b$Kmax*dd$dimb & prior.b$Kmax > init.b$K) init.b$Sigma <- init.b$Sigma[1:(init.b$K*dd$dimb),]
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


########## ========== tuneMCMC ========== ##########
########## ============================== ##########
  if (!is.list(tuneMCMC)) stop("tuneMCMC must be a list")
  intuneMCMC <- names(tuneMCMC)
  itune.beta <- match("beta", intuneMCMC, nomatch=NA)
  itune.b    <- match("b", intuneMCMC, nomatch=NA)
  
  if (dd$Rd){
    if (is.na(itune.beta)) tuneMCMC$beta <- rep(1, dd$Rd)
    if (length(tuneMCMC$beta) == 1) tuneMCMC$beta <- rep(tuneMCMC$beta, dd$Rd)
    if (length(tuneMCMC$beta) != dd$Rd) stop(paste("tuneMCMC$beta must be of length ", dd$Rd, sep=""))
    if (any(is.na(tuneMCMC$beta))) stop("NA in tuneMCMC$beta")        
    if (any(tuneMCMC$beta <= 0)) stop("tuneMCMC$beta must be all positive")            
  }else{
    tuneMCMC$beta <- 1
  }  
  
  if (dd$dimb){
    if (is.na(itune.b)) tuneMCMC$b <- 1
    if (length(tuneMCMC$b) != 1) stop(paste("tuneMCMC$b must be of length ", 1, sep=""))
    if (any(is.na(tuneMCMC$b))) stop("NA in tuneMCMC$b")        
    if (any(tuneMCMC$b <= 0)) stop("tuneMCMC$b must be all positive")                
  }else{
    tuneMCMC$b <- 1
  }  
  
  
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


########## ========== Some additional parameters ##########
########## ===================================== ##########
  if (prior.b$priorK == "fixed") lsum_Ir_b <- ifit$I * CK_b
  else                           lsum_Ir_b <- 1

  
####### Parameters passed to C++ which will be stored also in the resulting object (to be able to use them in related functions)
####### ========================================================================================================================
  R_cd <- c(dd$Rc, dd$Rd)
  names(R_cd) <- c("R_c", "R_d")
  p_fI_q_rI <- c(dd$p, dd$CfixedIntcpt, dd$q, dd$CrandomIntcpt)
  names(p_fI_q_rI) <- paste(rep(c("p", "fixedIntcpt", "q", "randomIntcpt"), each=sum(R_cd)), rep(1:sum(R_cd), 4), sep="")
  
  Cpar <- list(Y_c              = ifit$Cy_c,
               Y_d              = ifit$Cy_d,
               R_cd             = R_cd,
               dist             = dd$ndist,
               I                = ifit$I,
               n                = ifit$Cn,
               X                = ifit$CX,
               Z                = ifit$CZ,
               p_fI_q_rI        = p_fI_q_rI,
               priorDouble_eps  = peps$CpriorDouble_eps,
               priorInt_b       = pbb$CpriorInt_b,
               priorDouble_b    = pbb$CpriorDouble_b,
               priorDouble_beta = pbeta$CpriorDouble_beta,
               tune_scale_beta  = tuneMCMC$beta,
               tune_scale_b     = tuneMCMC$b)
  ifit$Cy_c <- NULL
  ifit$Cy_d <- NULL
  ifit$C_n  <- NULL
  ifit$CX   <- NULL
  ifit$CZ   <- NULL  

  #return(list(Cpar=Cpar, scb=scb, nMCMC=nMCMC, store=store,
  #            Csigma_eps=Csigma_eps, CgammaInv_eps=CgammaInv_eps,
  #            CK_b=CK_b, Cw_b=Cw_b, Cmu_b=Cmu_b,
  #            dd=dd, prior.b=prior.b,
  #            CLi_b=CLi_b, CgammaInv_b=CgammaInv_b, Cr_b=Cr_b,
  #            ifit=ifit, Cbeta=Cbeta, Cbb=Cbb, lsum_Ir_b=lsum_Ir_b))
    
########## ========== MCMC simulation                              ========== ##########
########## ================================================================== ##########
  cat(paste("MCMC sampling started on ", date(), ".\n", sep=""))  
  MCMC <- .C("GLMM_MCMC",
             Y_c                       = as.double(Cpar$Y_c),
             Y_d                       = as.integer(Cpar$Y_d),
             keepChain_nMCMC_R_cd_dist = as.integer(c(store, nMCMC, Cpar$R_cd, Cpar$dist)),
             I_n                       = as.integer(c(Cpar$I, Cpar$n)),
             X                         = as.double(Cpar$X),
             #XtX                       = as.double(ifit$CXtX),               ### REMOVED ON 21/10/2009, XtX is computed directly in C++
             Z                         = as.double(Cpar$Z),
             #ZitZi                     = as.double(ifit$CZitZi),             ### REMOVED ON 20/10/2009, ZitZi is computed directly in C++
             p_fI_q_rI                 = as.integer(Cpar$p_fI_q_rI),
             shiftScale_b              = as.double(scb$CshiftScale_b),
             priorDouble_eps           = as.double(Cpar$priorDouble_eps),
             priorInt_b                = as.integer(Cpar$priorInt_b),
             priorDouble_b             = as.double(Cpar$priorDouble_b),
             priorDouble_beta          = as.double(Cpar$priorDouble_beta),
             tune_scale_beta           = as.double(Cpar$tune_scale_beta),
             tune_scale_b              = as.double(Cpar$tune_scale_b),             
             sigma_eps                 = as.double(Csigma_eps),
             gammaInv_eps              = as.double(CgammaInv_eps),
             K_b                       = as.integer(CK_b),
             w_b                       = as.double(Cw_b),
             mu_b                      = as.double(Cmu_b),
             Q_b                       = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             Sigma_b                   = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             Li_b                      = as.double(CLi_b),
             gammaInv_b                = as.double(CgammaInv_b),
             r_b                       = as.integer(Cr_b),
             r_b_first                 = integer(ifit$I),
             beta                      = as.double(Cbeta),
             b                         = as.double(Cbb),
             b_first                   = double(length(Cbb)),
             chsigma_eps               = double(ifelse(dd$Rc, dd$Rc * nMCMC["keep"], 1)),
             chgammaInv_eps            = double(ifelse(dd$Rc, dd$Rc * nMCMC["keep"], 1)),
             chK_b                     = integer(ifelse(dd$dimb, nMCMC["keep"], 1)),
             chw_b                     = double(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chmu_b                    = double(ifelse(dd$dimb, dd$dimb * prior.b$Kmax * nMCMC["keep"], 1)),
             chQ_b                     = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chSigma_b                 = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chLi_b                    = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax * nMCMC["keep"], 1)),
             chgammaInv_b              = double(ifelse(dd$dimb, dd$dimb * nMCMC["keep"], 1)),
             chorder_b                 = integer(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chrank_b                  = integer(ifelse(dd$dimb, prior.b$Kmax * nMCMC["keep"], 1)),
             chMeanData_b              = double(ifelse(dd$dimb, dd$dimb * nMCMC["keep"], 1)),
             chCorrData_b              = double(ifelse(dd$dimb, dd$LTb * nMCMC["keep"], 1)),
             chbeta                    = double(ifelse(dd$lbeta, dd$lbeta * nMCMC["keep"], 1)),
             chb                       = double(ifelse(dd$dimb, ifelse(store["b"], ifit$I * dd$dimb * nMCMC["keep"], ifit$I * dd$dimb), 1)),
             chGLMMLogL                = double(nMCMC["keep"]),
             chLogL                    = double(nMCMC["keep"]),
             naccept_beta              = integer(dd$Rc + dd$Rd),
             naccept_b                 = integer(ifit$I),
             pm_eta_fixed              = double(ifit$sumCn),
             pm_eta_random             = double(ifit$sumCn),
             pm_meanY                  = double(ifit$sumCn),
             pm_stres                  = double(ifit$sumCn),
             pm_b                      = double(ifelse(dd$dimb, dd$dimb * ifit$I, 1)),
             pm_w_b                    = double(ifelse(dd$dimb, prior.b$Kmax, 1)),
             pm_mu_b                   = double(ifelse(dd$dimb, dd$dimb * prior.b$Kmax, 1)),
             pm_Q_b                    = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_Sigma_b                = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_Li_b                   = double(ifelse(dd$dimb, dd$LTb * prior.b$Kmax, 1)),
             pm_indGLMMLogL            = double(ifit$I),
             pm_indLogL                = double(ifit$I),
             pm_indLogpb               = double(ifit$I),
             sum_Ir_b                  = integer(lsum_Ir_b),
             sum_Pr_b_b                = double(lsum_Ir_b),
             iter                      = as.integer(0),
             err                       = as.integer(0),
             PACKAGE=thispackage)            
  cat(paste("MCMC sampling finished on ", date(), ".\n", sep=""))
  if (MCMC$err) stop("Something went wrong.")


  ########## ========== State of MCMC (last and first kept) ========== ##########
  ########## ========================================================= ##########
  if (dd$dimb){
    state.w_b       <- as.numeric(MCMC$w_b[1:MCMC$K_b])
    state_first.w_b <- as.numeric(MCMC$chw_b[1:MCMC$chK_b[1]])    
    names(state.w_b)       <- paste("w", 1:MCMC$K_b, sep="")
    names(state_first.w_b) <- paste("w", 1:MCMC$chK_b[1], sep="")    
    
    state.r_b       <- as.numeric(MCMC$r_b + 1)
    state_first.r_b <- as.numeric(MCMC$r_b_first + 1)    
    names(state.r_b) <- names(state_first.r_b) <- paste("r", 1:ifit$I, sep="")
    
    state.gammaInv_b       <- as.numeric(MCMC$gammaInv_b)
    state_first.gammaInv_b <- as.numeric(MCMC$chgammaInv_b[1:dd$dimb])    
    names(state_first.gammaInv_b) <- names(state.gammaInv_b) <- paste("gammaInv", 1:dd$dimb, sep="")    
    
    if (dd$dimb == 1){
      state.mu_b       <- as.numeric(MCMC$mu_b[1:MCMC$K_b])
      state_first.mu_b <- as.numeric(MCMC$chmu_b[1:MCMC$chK_b[1]])      
      names(state.mu_b)       <- paste("mu", 1:MCMC$K_b, sep="")
      names(state_first.mu_b) <- paste("mu", 1:MCMC$chK_b[1], sep="")      
      
      state.Li_b       <- as.numeric(MCMC$Li_b[1:MCMC$K_b])
      state_first.Li_b <- as.numeric(MCMC$chLi_b[1:MCMC$chK_b[1]])      
      names(state.Li_b)       <- paste("Li", 1:MCMC$K_b, sep="")
      names(state_first.Li_b) <- paste("Li", 1:MCMC$chK_b[1], sep="")      
      
      state.Sigma_b       <- (1 / state.Li_b)^2
      state_first.Sigma_b <- (1 / state_first.Li_b)^2      
      names(state.Sigma_b)       <- paste("Sigma", 1:MCMC$K_b, sep="")
      names(state_first.Sigma_b) <- paste("Sigma", 1:MCMC$chK_b[1], sep="")      
      
      state.Q_b       <- as.numeric(MCMC$Q_b[1:MCMC$K_b])
      state_first.Q_b <- as.numeric(MCMC$chQ_b[1:MCMC$chK_b[1]])      
      names(state.Q_b)       <- paste("Q", 1:MCMC$K_b, sep="")
      names(state_first.Q_b) <- paste("Q", 1:MCMC$chK_b[1], sep="")      

      state.b       <- as.numeric(MCMC$b)
      state_first.b <- as.numeric(MCMC$b_first)      
      names(state.b) <- names(state_first.b) <- 1:ifit$I
    }else{
      state.mu_b       <- matrix(MCMC$mu_b[1:(dd$dimb*MCMC$K_b)], ncol=dd$dimb, byrow=TRUE)
      state_first.mu_b <- matrix(MCMC$chmu_b[1:(dd$dimb*MCMC$chK_b[1])], ncol=dd$dimb, byrow=TRUE)      
      rownames(state.mu_b)       <- paste("j", 1:MCMC$K_b, sep="")
      rownames(state_first.mu_b) <- paste("j", 1:MCMC$chK_b[1], sep="")      
      colnames(state.mu_b) <- colnames(state_first.mu_b) <- paste("m", 1:dd$dimb, sep="")
      
      state.Li_b       <- as.numeric(MCMC$Li_b[1:(dd$LTb*MCMC$K_b)])
      state_first.Li_b <- as.numeric(MCMC$chLi_b[1:(dd$LTb*MCMC$chK_b[1])])      
      names(state.Li_b)       <- paste("Li", rep(1:MCMC$K_b, each=dd$LTb), rep(dd$naamLTb, MCMC$K_b), sep="")
      names(state_first.Li_b) <- paste("Li", rep(1:MCMC$chK_b[1], each=dd$LTb), rep(dd$naamLTb, MCMC$chK_b[1]), sep="")      
      
      state.Sigma_b       <- matrix(NA, ncol=dd$dimb, nrow=dd$dimb*MCMC$K_b)
      rownames(state.Sigma_b) <- paste("j", rep(1:MCMC$K_b, each=dd$dimb), ".", rep(1:dd$dimb, MCMC$K_b), sep="")
      colnames(state.Sigma_b) <- paste("m", 1:dd$dimb, sep="")      
      for (j in 1:MCMC$K_b){
        tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state.Li_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpSigma <- tmpSigma %*% t(tmpSigma)
        tmpSigma <- chol2inv(chol(tmpSigma))
        state.Sigma_b[((j-1)*dd$dimb+1):(j*dd$dimb),] <- tmpSigma
      }

      state_first.Sigma_b <- matrix(NA, ncol=dd$dimb, nrow=dd$dimb*MCMC$chK_b[1])      
      rownames(state_first.Sigma_b) <- paste("j", rep(1:MCMC$chK_b[1], each=dd$dimb), ".", rep(1:dd$dimb, MCMC$chK_b[1]), sep="")      
      colnames(state_first.Sigma_b) <- paste("m", 1:dd$dimb, sep="")      
      for (j in 1:MCMC$chK_b[1]){
        tmpSigma <- matrix(0, nrow=dd$dimb, ncol=dd$dimb)
        tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- state_first.Li_b[((j-1)*dd$LTb+1):(j*dd$LTb)]
        tmpSigma <- tmpSigma %*% t(tmpSigma)
        tmpSigma <- chol2inv(chol(tmpSigma))
        state_first.Sigma_b[((j-1)*dd$dimb+1):(j*dd$dimb),] <- tmpSigma
      }
      
      state.Q_b       <- as.numeric(MCMC$Q_b[1:(dd$LTb*MCMC$K_b)])
      state_first.Q_b <- as.numeric(MCMC$chQ_b[1:(dd$LTb*MCMC$chK_b[1])])      
      names(state.Q_b)       <- paste("Q", rep(1:MCMC$K_b, each=dd$LTb), rep(dd$naamLTb, MCMC$K_b), sep="")
      names(state_first.Q_b) <- paste("Q", rep(1:MCMC$chK_b[1], each=dd$LTb), rep(dd$naamLTb, MCMC$chK_b[1]), sep="")
      
      state.b       <- matrix(MCMC$b, ncol=dd$dimb, nrow=ifit$I, byrow=TRUE)
      state_first.b <- matrix(MCMC$b_first, ncol=dd$dimb, nrow=ifit$I, byrow=TRUE)      
      colnames(state.b) <- colnames(state_first.b) <- paste("b", 1:dd$dimb, sep="")
      rownames(state.b) <- rownames(state_first.b) <- 1:ifit$I
      
    }
    nCompTotal_b<- sum(MCMC$chK_b)
    freqK_b <- table(MCMC$chK_b)
    propK_b <- prop.table(freqK_b)
  }else{
    state.w_b <- state.r_b <- state.gamma_b <- state.mu_b <- state.Li_b <- state.Sigma_b <- state.Q_b <- state.b <- 0
    state_first.w_b <- state_first.r_b <- state_first.gamma_b <- state_first.mu_b <- state_first.Li_b <- state_first.Sigma_b <- state_first.Q_b <- state_first.b <- 0    
  }  

  if (dd$lbeta){
    state.beta       <- as.numeric(MCMC$beta)
    state_first.beta <- as.numeric(MCMC$chbeta[1:dd$lbeta])    
    names(state.beta) <- names(state_first.beta) <- paste("beta", 1:dd$lbeta, sep="")
  }else{
    state.beta <- state_first.beta <- 0
  }  

  if (dd$Rc){
    state.sigma_eps       <- as.numeric(MCMC$sigma_eps)
    state_first.sigma_eps <- as.numeric(MCMC$chsigma_eps[1:dd$Rc])    
    names(state.sigma_eps) <- names(state_first.sigma_eps) <- paste("sigma", 1:dd$Rc, sep="")
    
    state.gammaInv_eps <- as.numeric(MCMC$gammaInv_eps)
    state_first.gammaInv_eps <- as.numeric(MCMC$chgammaInv_eps[1:dd$Rc])    
    names(state.gammaInv_eps) <- names(state_first.gammaInv_eps) <- paste("gammaInv", 1:dd$Rc, sep="")
  }else{
    state.sigma_eps <- state.gammaInv_eps <- 0
    state_first.sigma_eps <- state_first.gammaInv_eps <- 0    
  }  

  
  ########## ========== Performance of MCMC ========== ##########
  ########## ========================================= ##########
  prop.accept.beta <- MCMC$naccept_beta / (nMCMC["keep"] * nMCMC["thin"])
  if (dd$R > 1) names(prop.accept.beta) <- NAAM.RESP
  prop.accept.b <- MCMC$naccept_b / (nMCMC["keep"] * nMCMC["thin"])  
    
  
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
              init.eps         = init.eps)


  if (dd$lbeta){
    RET$init.beta        <- init.beta
    RET$state.first.beta <- state_first.beta    
    RET$state.last.beta  <- state.beta
    RET$prop.accept.beta <- prop.accept.beta    
  }  
  
  if (dd$dimb){
    RET$init.b  <- init.b
    RET$state.first.b <- list(b        = state_first.b,               
                              K        = as.numeric(MCMC$chK_b[1]),  
                              w        = state_first.w_b,             
                              mu       = state_first.mu_b,            
                              Sigma    = state_first.Sigma_b,         
                              Li       = state_first.Li_b,            
                              Q        = state_first.Q_b,             
                              gammaInv = state_first.gammaInv_b,      
                              r        = state_first.r_b)    
    RET$state.last.b <- list(b        = state.b,               
                             K        = as.numeric(MCMC$K_b),  
                             w        = state.w_b,             
                             mu       = state.mu_b,            
                             Sigma    = state.Sigma_b,         
                             Li       = state.Li_b,            
                             Q        = state.Q_b,             
                             gammaInv = state.gammaInv_b,      
                             r        = state.r_b)
    RET$prop.accept.b <- prop.accept.b
    RET$scale.b <- scale.b
    RET$freqK_b <- freqK_b
    RET$propK_b <- propK_b    
  }                                 

  if (dd$Rc){
    RET$init.eps <- init.eps
    RET$state.first.eps <- list(sigma    = state_first.sigma_eps,
                                gammaInv = state_first.gammaInv_eps)
    RET$state.last.eps <- list(sigma    = state.sigma_eps,
                               gammaInv = state.gammaInv_eps)    
  }  
                                      

  ########## ========== Posterior means of quantities computed in C++ ========== ##########
  ########## =================================================================== ##########
  RET$poster.mean.y <- list()
  used <- 0
  s    <- 1
  while (s <= Cpar$R_cd["R_c"]){
    ns    <- Cpar$n[((s-1)*Cpar$I+1):(s*Cpar$I)]
    index <- (used+1):(used + sum(ns))    
    used  <- index[length(index)]
    RET$poster.mean.y[[s]] <- data.frame(id         = rep(1:Cpar$I, ns),
                                         observed   = Cpar$Y_c[index],
                                         fitted     = as.numeric(MCMC$pm_meanY[index]),
                                         stres      = as.numeric(MCMC$pm_stres[index]),
                                         eta.fixed  = as.numeric(MCMC$pm_eta_fixed[index]),
                                         eta.random = as.numeric(MCMC$pm_eta_random[index]))
    s <- s + 1
  }
  used2 <- 0
  while (s <= Cpar$R_cd["R_c"] + Cpar$R_cd["R_d"]){
    ns     <- Cpar$n[((s-1)*Cpar$I+1):(s*Cpar$I)]
    index  <- (used+1):(used + sum(ns))    
    used   <- index[length(index)]
    index2 <- (used2+1):(used2 + sum(ns))    
    used2  <- index2[length(index2)]  
    RET$poster.mean.y[[s]] <- data.frame(id         = rep(1:Cpar$I, ns),
                                         observed   = Cpar$Y_d[index2],
                                         fitted     = as.numeric(MCMC$pm_meanY[index]),
                                         stres      = as.numeric(MCMC$pm_stres[index]),
                                         eta.fixed  = as.numeric(MCMC$pm_eta_fixed[index]),
                                         eta.random = as.numeric(MCMC$pm_eta_random[index]))
    s <- s + 1
  }  
  names(RET$poster.mean.y) <- colnames(dd$y)
  
  if (dd$dimb){
    MCMC$pm_b <- matrix(MCMC$pm_b, ncol=dd$dimb, byrow=TRUE)
    RET$poster.mean.profile <- as.data.frame(MCMC$pm_b)
    colnames(RET$poster.mean.profile) <- paste("b", 1:dd$dimb, sep="")

    RET$poster.mean.profile$Logpb         <- as.numeric(MCMC$pm_indLogpb)    
    RET$poster.mean.profile$Cond.Deviance <- as.numeric(-2 * MCMC$pm_indLogL)
    RET$poster.mean.profile$Deviance      <- as.numeric(-2 * MCMC$pm_indGLMMLogL)
    
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
        RET$poster.mean.Li_b[[j]] <- tmpLi      
      }
      names(RET$poster.mean.Q_b) <- names(RET$poster.mean.Sigma_b) <- names(RET$poster.mean.Li_b) <- paste("j", 1:prior.b$Kmax, sep="")    
    }      
  }else{
    RET$poster.mean.profile <- data.frame(LogL     = as.numeric(MCMC$pm_indLogL),
                                          Deviance = as.numeric(-2 * MCMC$pm_indGLMMLogL))
  }  


  ########## ========== Clustering based on posterior P(alloc = k | y) or on P(alloc = k | theta, b, y)    ========== ##########
  ########## ======================================================================================================== ##########
  if (dd$dimb){
    if (prior.b$priorK == "fixed"){
      if (CK_b == 1){
        RET$poster.comp.prob1 <- RET$poster.comp.prob2 <- matrix(1, nrow = ifit$I, ncol = 1)
      }else{

        ### Using mean(I(r=k))
        MCMC$sum_Ir_b <- matrix(MCMC$sum_Ir_b, ncol = CK_b, nrow = ifit$I, byrow = TRUE)
        Denom <- apply(MCMC$sum_Ir_b, 1, sum)
        RET$poster.comp.prob1 <- MCMC$sum_Ir_b / matrix(rep(Denom, CK_b), ncol = CK_b, nrow = ifit$I)

        ### Using mean(P(r=k | theta, b, y))
        MCMC$sum_Pr_b_b<- matrix(MCMC$sum_Pr_b_b, ncol = CK_b, nrow = ifit$I, byrow = TRUE)
        RET$poster.comp.prob2 <- MCMC$sum_Pr_b_b/ matrix(rep(Denom, CK_b), ncol = CK_b, nrow = ifit$I)        
      }  
    }  
  }  
  
    
  ########## ========== Additional posterior summaries                ========== ##########
  ########## =================================================================== ##########
  qProbs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  nSumm <- c("Mean", "Std.Dev.", "Min.", "2.5%", "1st Qu.", "Median", "3rd Qu.", "97.5%", "Max.")

  mean.Deviance  <- -2 * mean(MCMC$chGLMMLogL, na.rm=TRUE)
  quant.Deviance <- 2 * quantile(-MCMC$chGLMMLogL, prob=qProbs, na.rm=TRUE)  
  sd.Deviance    <- 2 * sd(MCMC$chGLMMLogL, na.rm=TRUE)
  summ.Deviance  <-  c(mean.Deviance, sd.Deviance, quant.Deviance)
  mean.Cond.Deviance  <- -2 * mean(MCMC$chLogL, na.rm=TRUE)
  quant.Cond.Deviance <- 2 * quantile(-MCMC$chLogL, prob=qProbs, na.rm=TRUE)  
  sd.Cond.Deviance    <- 2 * sd(MCMC$chLogL, na.rm=TRUE)
  summ.Cond.Deviance  <-  c(mean.Cond.Deviance, sd.Cond.Deviance, quant.Cond.Deviance)  
  RET$summ.Deviance <- data.frame(Deviance = summ.Deviance, Cond.Deviance = summ.Cond.Deviance)
  rownames(RET$summ.Deviance) <- nSumm
    
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
      RET$summ.beta <- rbind(mean.beta, sd.beta, quant.beta)
      RET$summ.beta <- as.data.frame(RET$summ.beta) 
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
    RET$Deviance      <- as.numeric(-2 * MCMC$chGLMMLogL)
    RET$Cond.Deviance <- as.numeric(-2 * MCMC$chLogL)    
    
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
        RET$b <- matrix(MCMC$chb, ncol=dd$dimb*ifit$I, byrow=TRUE)
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

  
  ########## ========== Additional objects (added on 08/26/2010) ========== ##########
  ########## ============================================================== ##########
  RET$relabel_b <- list(type="mean", par=1)       #### default re-labeling is performed using the first margin of the mixture means
  RET$Cpar <- Cpar
  
  class(RET) <- "GLMM_MCMC"
  return(RET)    
}

