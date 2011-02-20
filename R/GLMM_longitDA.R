##
##  PURPOSE:   Longitudinal discriminant analysis
##             based on GLMM MCMC fits (with possibly several response variables) of 
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    05/08/2009 (as GLMM_longitClust)
##              28/10/2009: renamed to GLMM_longitDA
##
##  FUNCTIONS:  GLMM_longitDA
##
## ================================================================================================

## *************************************************************
## GLMM_longitDA
## *************************************************************
##
GLMM_longitDA <- function(mod, w.prior, y, id, time, x, z, xz.common=TRUE, info)
{
  thispackage <- "mixAK"

########## ========== Some input checks ========== ##########
########## ======================================= ##########
  if (!is.list(mod)) stop("mod must be a list")
  if (length(mod) < 2) stop("mod must be a list of length 2 or more")
  if (any(sapply(mod, class) != "GLMM_MCMC")) stop("all components of mod must be of class GLMM_MCMC")
  nClust <- length(mod)

  if (length(w.prior) != nClust) stop(paste("w.prior must be of length ", nClust, sep=""))
  if (any(w.prior < 0)) stop("w.prior must all be non-negative")
  w.prior <- w.prior / sum(w.prior)
  
  Rc <- mod[[1]]$R["Rc"]
  Rd <- mod[[1]]$R["Rd"]
  R <- Rc + Rd

  keepMCMC <- mod[[1]]$nMCMC["keep"]
  for (cl in 2:nClust){
    if (mod[[cl]]$R["Rc"] != Rc) stop(paste("cluster number ", cl, " has ", mod[[cl]]$R["Rc"], " cont. responses (must be ", Rc, ")", sep=""))
    if (mod[[cl]]$R["Rd"] != Rd) stop(paste("cluster number ", cl, " has ", mod[[cl]]$R["Rd"], " discrete responses (must be ", Rd, ")", sep=""))
    keepMCMC <- c(keepMCMC, mod[[cl]]$nMCMC["keep"])
  }  

  
########## ========== Design matrices ============ ##########
########## ======================================= ##########
  if (!xz.common){
    if (!is.list(x)) stop("x must be a list when xz.common is FALSE")
    if (length(x) != nClust) stop(paste("x must be a list of length ", nClust, sep=""))
    if (!is.list(z)) stop("z must be a list when xz.common is FALSE")
    if (length(z) != nClust) stop(paste("z must be a list of length ", nClust, sep=""))    
  }  


########## ========== Work out the data ============ ##########
########## ========================================= ##########
  Cp <- Cq <- CfixedIntcpt <- CrandomIntcpt <- numeric()
  CX <- CZ <- numeric()
  ## CXtX <- CZitZi <- numeric()                              ### CODE REMOVED ON 02/11/2009
  CshiftScale_b <- numeric()
  Kmax_b <- numeric()
  chK_b <- chw_b <- chmu_b <- chLi_b <- numeric()
  chsigma_eps <- numeric()
  chalpha <- numeric()
  
  for (cl in 1:nClust){
    if (xz.common){
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x, z=z, random.intercept=mod[[cl]]$random.intercept)
    }else{
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x[[cl]], z=z[[cl]], random.intercept=mod[[cl]]$random.intercept)
    }

    ifit <- GLMM_MCMCifit(do.init=FALSE, na.complete=TRUE,
                          y=dd$y, dist=dd$dist, id=dd$id, time=dd$time, x=dd$x, z=dd$z, random.intercept=dd$random.intercept,
                          xempty=dd$xempty, zempty=dd$zempty, Rc=dd$Rc, Rd=dd$Rd,
                          p=dd$p, p_fi=dd$p_fi, q=dd$q, q_ri=dd$q_ri, lalpha=dd$lalpha, dimb=dd$dimb)
    ## ifit$x[[s]], ifit$z[[s]] now contains also intercept columns if these should be included
    ## * also rows corresponding to missing values are removed
    ## * there is the same number of observations for each cluster
    
    ##### Variables common to all clusters
    if (cl == 1){
      dist1 <- dd$dist
      Cdist <- dd$ndist

      Cy_c  <- ifit$Cy_c
      Cy_d  <- ifit$Cy_d
      I     <- ifit$I
      Cn    <- ifit$Cn
      sumCn <- ifit$sumCn

      ID <- ifit$ID[[1]]
      TIME <- ifit$time

      nrow_pi <- sum(ifit$n[[1]])
    }  
    
    ##### Check whether supplied y is compatible with these used to fit mod
    if (ncol(dd$y) != R) stop(paste("y must have ", R, " columns", sep=""))

    ##### Check whether common variables (across clusters) are really common
    if (any(dd$dist != dist1)) stop(paste("dist of cluster ", cl, " is incompatible with that of cluster 1", sep=""))

    if (length(ifit$Cy_c) != length(Cy_c)) stop(paste("number of NA's for cluster ", cl, " leads to incompatibility in Cy_c", sep=""))
    if (length(ifit$Cy_d) != length(Cy_d)) stop(paste("number of NA's for cluster ", cl, " leads to incompatibility in Cy_d", sep=""))    
    if (ifit$I != I) stop(paste("I for cluster ", cl, " is incompatible with that of clulster 1", sep=""))
    if (length(ifit$Cn) != length(Cn)) stop(paste("Cn for cluster ", cl, " is incompatible with that of cluster 1", sep=""))
    if (any(ifit$Cn != Cn)) stop(paste("Cn for cluster ", cl, " is incompatible with that of cluster 1", sep=""))    
    
    ##### Check whether supplied x and z is compatible with these used to fit mod
    ##### Double check that there is the same number of observations for each response
    for (s in 1:R){      
      if (dd$xempty[s]){
        if (mod[[cl]]$p[[s]] != 0) stop(paste("x[[", s, "]] should not be empty for cluster ", cl, sep=""))
      }else{
        if (ncol(dd$x[[s]]) != mod[[cl]]$p[s]) stop(paste("x[[", s, "]] matrix for cluster ", cl, " must have ", mod[[cl]]$p[s], " columns", sep=""))
      }

      if (dd$zempty[s]){
        if (dd$z[[s]] != "empty") stop(paste("z[[", s, "]] must be empty for cluster ", cl, sep=""))
      }else{
        if (ncol(dd$z[[s]]) != mod[[cl]]$q[s]) stop(paste("z[[", s, "]] matrix for cluster ", cl, " must have ", mod[[cl]]$q[s], " columns", sep=""))
      }

      if (any(ifit$n[[s]] != ifit$n[[1]])) stop("BUG in the function, contact AK")
    }

    ### ===== CODE REMOVED ON 02/11/2009 =====
    ##### Calculate CZitZi to be used here (it is different from that needed by GLMM_MCMC function)
    ##### and now included in ifit$CZitZi
    ##### REMARK:  Due to a loop, this is relatively slow for larger data, do it in C++ in a future???
    ##ZitZi <- numeric(0)
    ##cumn <- c(0, cumsum(ifit$n[[1]]))
    ##for (i in 1:I){                       ## loop over longitudinal profiles
    ##  for (j in 1:ifit$n[[1]][i]){        ## loop over observations within longitudinal profiles
    ##    for (s in 1:R){
    ##      if (dd$q_ri[s]){
    ##        tmpZ <- matrix(ifit$z[[s]][(cumn[i]+1):(cumn[i]+j),], ncol=ncol(ifit$z[[s]]))
    ##        tmpZtZ <- t(tmpZ) %*% tmpZ
    ##        ZitZi <- c(ZitZi, tmpZtZ[lower.tri(tmpZtZ, diag=TRUE)])
    ##      }  
    ##    }          
    ##  }  
    ##}
    ##cat("length ZitZi[", cl, "] = ", length(ZitZi), "\n", sep="")
    ### ===== END OF CODE REMOVED ON 02/11/2009 =====
    
    ##### Vectors to be supplied to C++
    Cp <- c(Cp, dd$p)
    Cq <- c(Cq, dd$q)    
    CfixedIntcpt  <- c(CfixedIntcpt, dd$CfixedIntcpt)
    CrandomIntcpt <- c(CrandomIntcpt, dd$CrandomIntcpt)
    CX     <- c(CX, ifit$CX)
    ##CXtX   <- c(CXtX, ifit$CXtX)           ## CODE REMOVED ON 02/11/2009
    CZ     <- c(CZ, ifit$CZ)
    ##CZitZi <- c(CZitZi, ZitZi)             ## CODE REMOVED ON 02/11/2009

    CshiftScale_b <- c(CshiftScale_b, mod[[cl]]$scale.b$shift, mod[[cl]]$scale.b$scale)
    
    Kmax_b <- c(Kmax_b, mod[[cl]]$prior.b$Kmax)
    if (mod[[cl]]$prior.b$priorK == "fixed"){
      chK_b  <- c(chK_b, mod[[cl]]$K_b)
      chw_b  <- c(chw_b, t(mod[[cl]]$w_b))
      chmu_b <- c(chmu_b, t(mod[[cl]]$mu_b))    
      chLi_b <- c(chLi_b, t(mod[[cl]]$Li_b))
    }else{
      stop("not implemented for variable K")
    }  
    
    if (mod[[cl]]$lalpha)   chalpha <- c(chalpha, t(mod[[cl]]$alpha))
    if (mod[[cl]]$R["Rc"]) chsigma_eps <- c(chsigma_eps, t(mod[[cl]]$sigma_eps))
  }

  if (!length(chalpha))      chalpha <- 0
  if (!length(chsigma_eps)) chsigma_eps <- 0  
  
  ### !!! CZ, CXtX, CZ, CZitZi contain one 0 when there is no such matrix
  ###     in a model for specific cluster

  if (missing(info)) info <- min(keepMCMC)
  if (info <= 0) info <- 1

  
########## ========== Validation of C++ code  ============ ##########
########## =============================================== ##########  
  DEBUG <- FALSE
  if (DEBUG){         ### Compute prediction in one iteration for selected observation completely in R
    cl <- 2
    iter <- 67
    ii <- 2
    jj <- 9
    cat("\nShowing R expressions for cl = ", cl, ", iter = ", iter, ", i = ", ii, ", j = ", jj, "\n================================================================\n", sep="")
    
    if (xz.common){
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x, z=z, random.intercept=mod[[cl]]$random.intercept)
    }else{
      dd <- GLMM_MCMCdata(y=y, dist=mod[[cl]]$dist, id=id, time=time, x=x[[cl]], z=z[[cl]], random.intercept=mod[[cl]]$random.intercept)
    }
    ifit <- GLMM_MCMCifit(do.init=FALSE, na.complete=TRUE,
                          y=dd$y, dist=dd$dist, id=dd$id, time=dd$time, x=dd$x, z=dd$z, random.intercept=dd$random.intercept,
                          xempty=dd$xempty, zempty=dd$zempty, Rc=dd$Rc, Rd=dd$Rd,
                          p=dd$p, p_fi=dd$p_fi, q=dd$q, q_ri=dd$q_ri, lalpha=dd$lalpha, dimb=dd$dimb)    
    
    LLTb <- (mod[[cl]]$dimb * (mod[[cl]]$dim + 1)) / 2
    
    ss <- mod[[cl]]$scale.b$shift
    SS <- diag(mod[[cl]]$scale.b$scale)

    balpha <- mod[[cl]]$alpha[iter,]
    
    KK  <- mod[[cl]]$K_b[iter]
    ww  <- mod[[cl]]$w_b[iter,]
    mmu <- matrix(mod[[cl]]$mu_b[iter,], nrow=KK, byrow=TRUE)
    LLi <- iLLi <- QQ <- DD <- list()
    QQmmu <- matrix(0, nrow=KK, ncol=mod[[cl]]$dimb)
    llog_w_EB <- numeric(KK) 
    for (k in 1:KK){
      LLi[[k]] <- diag(mod[[cl]]$dimb)
      LLi[[k]][lower.tri(LLi[[k]], diag=TRUE)] <- mod[[cl]]$Li_b[iter, ((k-1)*LLTb + 1):(k*LLTb)]
      iLLi[[k]] <- solve(LLi[[k]])
      QQ[[k]] <- LLi[[k]] %*% t(LLi[[k]])
      DD[[k]] <- chol2inv(t(LLi[[k]]))
      QQmmu[k,] <- as.numeric(QQ[[k]] %*% mmu[k,])
      llog_w_EB[k] <- log(ww[k]) - 0.5 * log(det(DD[[k]])) - 0.5 * t(mmu[k,]) %*% QQ[[k]] %*% mmu[k,]
    }
    names(llog_w_EB) <- rownames(QQmmu) <- rownames(mmu) <- names(LLi) <- names(iLLi) <- names(QQ) <- names(DD) <- paste("k=", 1:KK, sep="")
    colnames(QQmmu) <- colnames(mmu) <- 1:mod[[cl]]$dimb
    
    iID <- ifit$ID[[1]] == unique(ifit$ID[[1]])[ii]
    YYi <- ifit$Y[[1]][iID][1:jj]
    if (R > 1){
      for (s in 2:R) YYi <- c(YYi, ifit$Y[[s]][iID][1:jj])
    }
    SSigma <- diag(rep(mod[[cl]]$sigma_eps[iter,]^2, each=jj))
    iSSigma <- diag(rep(1 / mod[[cl]]$sigma_eps[iter,]^2, each=jj))    

    ZZ <- matrix(ifit$z[[1]][iID,], ncol=ncol(ifit$z[[1]]))
    if (R > 1){
      ZZ <- rbind(matrix(ZZ[1:jj,], ncol=ncol(ifit$z[[1]])), matrix(0, nrow=jj*(R-1), ncol=ncol(ZZ)))
      for (s in 2:R){
        ZZs <- matrix(ifit$z[[s]][iID,], ncol=ncol(ifit$z[[s]]))
        if (s == R){
          ZZs <- rbind(matrix(0, nrow=jj*(R-1), ncol=ncol(ZZs)), matrix(ZZs[1:jj,], ncol=ncol(ifit$z[[s]])))
        }else{
          ZZs <- rbind(matrix(0, nrow=jj*(s-1), ncol=ncol(ZZs)), matrix(ZZs[1:jj,], ncol=ncol(ifit$z[[s]])), matrix(0, nrow=jj*(R-s), ncol=ncol(ZZs)))
        }  
        ZZ <- cbind(ZZ, ZZs)
      }  
    }else{
      ZZ <- ZZ[1:jj,]
    }

    if (dd$xempty[1]) stop("DEBUG code not implemented when there is no X matrix for response 1")      
    XX <- matrix(ifit$x[[1]][iID,], ncol=ncol(ifit$x[[1]]))
    if (R > 1){
      XX <- rbind(matrix(XX[1:jj,], ncol=ncol(ifit$x[[1]])), matrix(0, nrow=jj*(R-1), ncol=ncol(XX)))
      for (s in 2:R){
        if (!dd$xempty[s]){
          XXs <- matrix(ifit$x[[s]][iID,], ncol=ncol(ifit$x[[s]]))
          if (s == R){
            XXs <- rbind(matrix(0, nrow=jj*(R-1), ncol=ncol(XXs)), matrix(XXs[1:jj,], ncol=ncol(ifit$x[[s]])))
          }else{
            XXs <- rbind(matrix(0, nrow=jj*(s-1), ncol=ncol(XXs)), matrix(XXs[1:jj,], ncol=ncol(ifit$x[[s]])), matrix(0, nrow=jj*(R-s), ncol=ncol(XXs)))
          }  
          XX <- cbind(XX, XXs)
        }
      }        
    }else{
      XX <- XX[1:jj,]
    }  
    
    ### Precisions, canonical means and means of full conditional distributions    
    AA <- list()
    cc <- bbhat <- bbhat2 <- matrix(0, nrow=KK, ncol=mod[[cl]]$dimb)
    mmu_full_part <- as.numeric(t(SS) %*% t(ZZ) %*% iSSigma %*% (YYi - XX %*% balpha - ZZ %*% ss))

    llog_w_EB_ij <- numeric(KK)
    for (k in 1:KK){
      tmp <- t(SS) %*% t(ZZ) %*% iSSigma
      AA[[k]] <- tmp %*% ZZ %*% SS + QQ[[k]]
      cc[k,] <- as.numeric(tmp %*% (YYi - XX %*% balpha - ZZ %*% ss) + QQ[[k]] %*% mmu[k,])
      bbhat[k,] <- solve(AA[[k]], cc[k,])
      chAA <- chol(AA[[k]])      
      bbhat2[k,] <- forwardsolve(t(chAA), cc[k,])
      bbhat2[k,] <- backsolve(chAA, bbhat2[k,])
      llog_w_EB_ij[k] <- llog_w_EB[k] + 0.5 * crossprod(cc[k,], bbhat[k,])
    }
    names(AA) <- rownames(cc) <- rownames(bbhat) <- rownames(bbhat2) <- names(llog_w_EB_ij) <- paste("k=", 1:KK, sep="")
    colnames(cc) <- colnames(bbhat) <- colnames(bbhat2) <- 1:mod[[cl]]$dimb    

    ww_EB_ij <- exp(llog_w_EB_ij - max(llog_w_EB_ij))
    ww_EB_ij <- ww_EB_ij / sum(ww_EB_ij)
    EEBscaled <- apply(ww_EB_ij * bbhat, 2, sum)
    EEB <- as.numeric(ss + SS %*% EEBscaled)
    EEBscaled2 <- apply(ww_EB_ij * bbhat2, 2, sum)
    EEB2 <- as.numeric(ss + SS %*% EEBscaled2)
    names(EEB) <- names(EEB2) <- 1:mod[[cl]]$dimb

    VV <- list()
    mm_marg <- matrix(0, nrow=KK, ncol=length(YYi))
    ff_b <- ff_y <- 0
    for (k in 1:KK){
      ff_b <- ff_b + ww[k] * dMVN(EEBscaled, mean=mmu[k,], Sigma=DD[[k]])
      VV[[k]] <- ZZ %*% SS %*% DD[[k]] %*% t(SS) %*% t(ZZ) + SSigma
      mm_marg[k,] <- as.numeric(XX %*% balpha + ZZ %*% ss + ZZ %*% SS %*% mmu[k,])
      ff_y <- ff_y + ww[k] * dMVN(YYi, mean=mm_marg[k,], Sigma=VV[[k]])      
    }
    names(VV) <- rownames(mm_marg) <- paste("k=", 1:KK, sep="")
    Ey_b <- as.numeric(XX %*% balpha + ZZ %*% EEB)
    ff_y_b <- dMVN(YYi, mean=Ey_b, Sigma=SSigma)
    
    #cat("Mixture precision matrices (Q matrices):\n")
    #print(QQ)
    #cat("Mixture matrices Li^{-1}:\n")
    #print(iLLi)    
    #cat("Mixture vectors Q %*% mu:\n")
    #print(QQmmu)
    #cat("Mixture log_w_EB (common part):\n")
    #print(llog_w_EB)

    ##cat("mu_full_part:\n")
    ##print(mmu_full_part)
    
    #cat("Precisions of full conditional distributions (A matrices):\n")
    #print(AA)
    #cat("Canonical means of full conditional distributions (c vectors):\n")
    #print(cc)
    #cat("\nMeans of full conditional distributions (bhat[k] vectors):\n")
    #print(bbhat)
    #print(bbhat2)    
    #cat("Mixture log_w_EB_ij (obs. specific):\n")
    #print(llog_w_EB_ij)
    #cat("Mixture w_EB_ij (obs. specific):\n")
    #print(ww_EB_ij)    
    #cat("EB(scaled):\n")
    #print(EEBscaled)
    #print(EEBscaled2)    
    #cat("EB:\n")
    #print(EEB)
    #print(EEB2)
    #cat("Marginal covariance matrices (V matrices):\n")
    #print(VV)
    #cat("Marginal means:\n")
    #print(mm_marg)        
    cat("f_b = ", ff_b, "\nf_y_b = ", ff_y_b, "\nf_y = ", ff_y, "\n\n", sep="")   
  }  

  if (FALSE){
    cat("Cy_c:\n")
    print(Cy_c)

    cat("Cy_d:\n")
    print(Cy_d)

    cat("CX:\n")
    print(CX)

    cat("CZ:\n")
    print(CZ)

    ##cat("CZitZi:\n")
    ##print(CZitZi)

    cat("R_c = ", Rc, ", R_d = ", Rd, ", nClust = ", nClust, ", I = ", I, "\n", sep="")
    cat("Cdist:\n")
    print(Cdist)

    cat("Cn:\n")
    print(Cn)

    cat("Cp:\n")
    print(Cp)
    
    cat("CfixedIntcpt:\n")
    print(CfixedIntcpt)

    cat("Cq:\n")
    print(Cq)
    
    cat("CrandomIntcpt:\n")
    print(CrandomIntcpt)

    cat("CshiftScale_b:\n")
    print(CshiftScale_b)
        
    cat("keepMCMC:\n")
    print(keepMCMC)        
  }  
  
  fit <- .C("GLMM_longitDA",
            Y_c          = as.double(Cy_c),
            R_c          = as.integer(Rc),
            Y_d          = as.integer(Cy_d),
            R_d          = as.integer(Rd),
            dist         = as.integer(Cdist),
            nClust       = as.integer(nClust),
            I            = as.integer(I),
            n            = as.integer(Cn),
            X            = as.double(CX),
            p            = as.integer(Cp),
            fixedIntcpt  = as.integer(CfixedIntcpt),
            Z            = as.double(CZ),
            #SZitZiS      = as.double(CZitZi),                  ### REMOVED ON 02/11/2009
            q            = as.integer(Cq),
            randIntcpt   = as.integer(CrandomIntcpt),
            shiftScale_b = as.double(CshiftScale_b),
            keepMCMC     = as.integer(keepMCMC),
            info         = as.integer(info),
            Kmax_b       = as.integer(Kmax_b),
            chsigma_eps  = as.double(chsigma_eps),
            chK_b        = as.integer(chK_b),
            chw_b        = as.double(chw_b),
            chmu_b       = as.double(chmu_b),
            chLi_b       = as.double(chLi_b),
            chalpha      = as.double(chalpha),
            pi_marg      = double(nrow_pi * nClust),
            pi_cond      = double(nrow_pi * nClust),
            pi_ranef     = double(nrow_pi * nClust),
            err          = as.integer(0),
            PACKAGE = thispackage)

  #browser()
  
  if (nrow_pi == 1){
    fit$pi_ranef <- w.prior * fit$pi_ranef
    fit$pi_ranef <- fit$pi_ranef / sum(fit$pi_ranef)

    fit$pi_cond <- w.prior * fit$pi_cond
    fit$pi_cond <- fit$pi_cond / sum(fit$pi_cond)

    fit$pi_marg <- w.prior * fit$pi_marg
    fit$pi_marg <- fit$pi_marg / sum(fit$pi_marg)    
  }else{     
    w.priorMat <- matrix(rep(w.prior, nrow_pi), ncol=nClust, byrow=TRUE)
    #
    fit$pi_ranef <- w.priorMat * matrix(fit$pi_ranef, ncol=nClust)
    fit$pi_ranef <- fit$pi_ranef / apply(fit$pi_ranef, 1, sum)
    #
    fit$pi_cond <- w.priorMat * matrix(fit$pi_cond, ncol=nClust)
    fit$pi_cond <- fit$pi_cond / apply(fit$pi_cond, 1, sum)
    #
    fit$pi_marg <- w.priorMat * matrix(fit$pi_marg, ncol=nClust)
    fit$pi_marg <- fit$pi_marg / apply(fit$pi_marg, 1, sum)
  }  

  RET <- list(ident = data.frame(id=ID, time=TIME),
              marg  = fit$pi_marg,
              cond  = fit$pi_cond,
              ranef = fit$pi_ranef)
  
  return(RET)  
}
