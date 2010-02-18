##
##  PURPOSE:   Re-labeling of the MCMC output.
##             * method for objects of class NMixMCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   08/02/2010
##
##  FUNCTION:  NMixRelabel.NMixMCMC (08/02/2010) 
##
## ======================================================================

## *************************************************************
## NMixRelabel.NMixMCMC
## *************************************************************
NMixRelabel.NMixMCMC <- function(object, type=c("mean", "weight", "stephens"), par, info, ...)
{
  thispackage <- "mixAK"

  LTp <- object$dim * (object$dim + 1)/2
  n <- object$Cpar$dimy["n"]
  
  ##### Determine re-labeling algorithm to use and additional parameters
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  type <- match.arg(type)
  if (type == "mean"){
    Ctype <- 1
    if (missing(par)) par <- 1
    if (!is.numeric(par)) stop("par must be a number")
    par <- par[1]
    iparam <- par
    if (iparam <= 0 | iparam > object$dim) stop(paste("par must be between 1 and ", object$dim, sep=""))

    object$relabel <- list(type=type, par=iparam)    ## resulting re-labeling
    iparam <- iparam - 1                             ## R -> C indexing
  }else{
    if (type == "weight"){
      Ctype <- 2
      par <- 1
      iparam <- 0
      object$relabel <- list(type=type, par=0)       ## resulting re-labeling
    }
    else{
      if (type == "stephens"){
        Ctype <- 3
        if (missing(par)) par <- list(type.init = "identity", par = 0, maxiter = 50)
        if (!is.list(par)) stop("par must be a list")
        inpar <- names(par)
        ind_type.init  <- match("type.init", inpar, nomatch=NA)
        ind_par        <- match("par",       inpar, nomatch=NA)
        ind_maxiter    <- match("maxiter",   inpar, nomatch=NA)

        if (is.na(ind_type.init)) par$type.init <- "identity"
        Ctype.init <- pmatch(par$type.init, table=c("identity", "mean", "weight"), nomatch=NA) - 1
        if (is.na(Ctype.init)) stop("unknown par$type.init option supplied")

        if (Ctype.init == 0){               ### par$init == identity
          par$par <- 1
        }else{
          if (Ctype.init == 1){             ### par$init == mean
            if (is.na(ind_par)) par$par <- 1
            if (!is.numeric(par$par)) stop("par$par must be a number")
            par$par <- par$par[1]
            if (par$par <= 0 | par$par > object$dim) stop(paste("par$par must be between 1 and ", object$dim, sep=""))            
          }else{
            if (Ctype.init == 2){           ### par$init == weight
              par$par <- 1
            }  
          }  
        }

        if (is.na(ind_maxiter)) par$maxiter <- 50
        if (!is.numeric(par$maxiter)) stop("par$maxiter must be a number")        
        par$maxiter <- par$maxiter[1]
        if (par$maxiter <= 0) stop("par$maxiter must be positive")

        iparam <- c(Ctype.init, par$par - 1, par$maxiter, 1)              ### iparam[4] determines whether search or transportation problem should
                                                                          ### be used in step2 (0 = transportation problem, 1 = search)
        names(iparam) <- c("type.init", "par", "maxiter", "type.step2")
      }  
    }
  }


  ##### Parameters of MCMC
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  if (is.null(object$K) | is.null(object$w) | is.null(object$mu) | is.null(object$Li) | is.null(object$Q) | is.null(object$Sigma)){
    stop("object does not contain sampled values")
  }  
  
  keepMCMC <- length(object$w) / object$K[1]
  if (missing(info)) info <- keepMCMC
  if (info <= 0 | info > keepMCMC) info <- keepMCMC

  
  ##### Some input checks
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  if (object$prior$priorK != "fixed") stop("only implemented for models with a fixed number of mixture components")


  ##### Initial values of censored observations
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  zinit <- (object$state.first$y - matrix(rep(object$scale$shift, n), ncol=object$dim, byrow=TRUE)) / matrix(rep(object$scale$scale, n), ncol=object$dim, byrow=TRUE)

  
  ##### Perform re-labeling
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  l_nchange <- ifelse(Ctype <= 2, 1, par$maxiter)
  
  MCMC <- .C("NMix_NMixRelabel",
             type         = as.integer(Ctype),
             iparam       = as.integer(iparam),
             y0           = as.double(t(object$Cpar$z0)),
             y1           = as.double(t(object$Cpar$z1)),
             censor       = as.integer(t(object$Cpar$censor)),
             dimy         = as.integer(object$Cpar$dimy),
             keepMCMC     = as.integer(keepMCMC),
             info         = as.integer(info),
             K            = as.integer(object$K[1]),
             chw          = as.double(t(object$w)),
             chmu         = as.double(t(object$mu)),
             chQ          = as.double(t(object$Q)),
             chSigma      = as.double(t(object$Sigma)),
             chLi         = as.double(t(object$Li)),
             chorder      = integer(object$K[1] * keepMCMC),
             chrank       = integer(object$K[1] * keepMCMC),
             y            = as.double(t(zinit)),
             r            = integer(n),
             pm_w         = double(object$K[1]),
             pm_mu        = double(object$dim * object$K[1]),
             pm_Q         = double(LTp * object$K[1]),
             pm_Sigma     = double(LTp * object$K[1]),
             pm_Li        = double(LTp * object$K[1]),
             sum_Ir       = integer(n * object$K[1]),
             hatPr_y      = double(n * object$K[1]),
             iter_relabel = as.integer(0),
             nchange      = integer(l_nchange),
             err          = as.integer(0),
             PACKAGE = thispackage)
  if (MCMC$err) stop("Something went wrong.")


  ##### New chains for order and rank (corresponding to newly labeled sample)
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$order <- matrix(as.numeric(MCMC$chorder + 1), ncol=object$K[1], byrow=TRUE)
  colnames(object$order) <- paste("order", 1:object$K[1], sep="")  
  MCMC$chorder <- NULL
  
  object$rank <- matrix(as.numeric(MCMC$chrank + 1), ncol=object$K[1], byrow=TRUE)
  colnames(object$rank) <- paste("rank", 1:object$K[1], sep="")                          
  MCMC$chrank <- NULL


  ##### Clustering based on posterior P(alloc = k | y) or on P(alloc = k | theta, b, y) 
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (object$K[1] == 1){
    object$poster.comp.prob1 <- object$poster.comp.prob2 <- matrix(1, nrow = n, ncol = 1)
  }else{

    ### Using mean(I(r=k))
    MCMC$sum_Ir <- matrix(MCMC$sum_Ir, ncol = object$K[1], nrow = n, byrow = TRUE)
    Denom <- apply(MCMC$sum_Ir, 1, sum)       ### this should be a vector of length n with all elements equal to the number of saved MCMC iterations 
    object$poster.comp.prob1 <- MCMC$sum_Ir / matrix(rep(Denom, object$K[1]), ncol = object$K[1], nrow = n)

    ### Using mean(P(r=k | theta, b, y))
    object$poster.comp.prob2 <- matrix(MCMC$hatPr_y, ncol = object$K[1], nrow = n, byrow = TRUE)
  }  

                        
  ##### Posterior means for mixture components
  ##### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$poster.mean.w <- as.numeric(MCMC$pm_w)
  names(object$poster.mean.w) <- paste("w", 1:object$K[1], sep="")

  object$poster.mean.mu <- matrix(MCMC$pm_mu, nrow=object$K[1], ncol=object$dim, byrow=TRUE)
  rownames(object$poster.mean.mu) <- paste("j", 1:object$K[1], sep="")
  colnames(object$poster.mean.mu) <- paste("m", 1:object$dim, sep="")

  object$poster.mean.Q <- object$poster.mean.Sigma <- object$poster.mean.Li <- list()
  for (j in 1:object$K[1]){
    tmpQ <- matrix(0, nrow=object$dim, ncol=object$dim)
    tmpQ[lower.tri(tmpQ, diag=TRUE)] <- MCMC$pm_Q[((j-1)*LTp+1):(j*LTp)]
    tmpQ[upper.tri(tmpQ, diag=FALSE)] <- t(tmpQ)[upper.tri(t(tmpQ), diag=FALSE)]
    object$poster.mean.Q[[j]] <- tmpQ
    
    tmpSigma <- matrix(0, nrow=object$dim, ncol=object$dim)
    tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- MCMC$pm_Sigma[((j-1)*LTp+1):(j*LTp)]
    tmpSigma[upper.tri(tmpSigma, diag=FALSE)] <- t(tmpSigma)[upper.tri(t(tmpSigma), diag=FALSE)]
    object$poster.mean.Sigma[[j]] <- tmpSigma
    
    tmpLi <- matrix(0, nrow=object$dim, ncol=object$dim)
    tmpLi[lower.tri(tmpLi, diag=TRUE)] <- MCMC$pm_Li[((j-1)*LTp+1):(j*LTp)]
    tmpLi[upper.tri(tmpLi, diag=FALSE)] <- t(tmpLi)[upper.tri(t(tmpLi), diag=FALSE)]
    object$poster.mean.Li[[j]] <- tmpLi      
  }
  names(object$poster.mean.Q) <- names(object$poster.mean.Sigma) <- names(object$poster.mean.Li) <- paste("j", 1:object$K[1], sep="")    


  ##### Numbers of changes of labelling at each re-labelling iteration (for Stephens' algorithm)
  ##### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #if (type == "stephens"){
  #  nchange <- MCMC$nchange[1:MCMC$iter_relabel]
  #  cat("Numbers of changes of labelling at each re-labelling iteration:\n")
  #  print(nchange)
  #}
  
  return(object)
}