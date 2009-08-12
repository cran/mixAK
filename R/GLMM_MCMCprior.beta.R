##
##  PURPOSE:   Handle prior.beta argument of GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    05/08/2009
##
##  FUNCTIONS:  GLMM_MCMCprior.beta
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCprior.beta
## *************************************************************
##
GLMM_MCMCprior.beta <- function(prior.beta, lbeta)
{
##### Variables in the resulting object:
#####          prior.beta  
#####          CpriorDouble_beta
##### -----------------------------------------------------------------------------------------------------------    
  if (lbeta){
    if (missing(prior.beta)) prior.beta <- list()
    if (!is.list(prior.beta)) stop("prior.beta must be a list")

    inprior.beta <- names(prior.beta)
    ibeta.mean   <- match("mean", inprior.beta, nomatch=NA)
    ibeta.var    <- match("var", inprior.beta, nomatch=NA)    

    ##### prior.beta:  mean
    ##### -----------------------------------------------
    if (is.na(ibeta.mean)) prior.beta$mean <- rep(0, lbeta)
    if (length(prior.beta$mean) == 1) prior.beta$mean <- rep(prior.beta$mean, lbeta)
    if (length(prior.beta$mean) != lbeta) stop(paste("prior.beta$mean must be of length", lbeta))
    if (any(is.na(prior.beta$mean))) stop("NA in prior.beta$mean")        
    Cbetamean <- as.numeric(prior.beta$mean)
    names(Cbetamean) <- names(prior.beta$mean) <- paste("beta", 1:lbeta, ".mean", sep="")
    
    ##### prior.beta:  var
    ##### -----------------------------------------------
    if (is.na(ibeta.var)) prior.beta$var <- rep(10000, lbeta)
    if (length(prior.beta$var) == 1) prior.beta$var <- rep(prior.beta$var, lbeta)
    if (length(prior.beta$var) != lbeta) stop(paste("prior.beta$var must be of length", lbeta))
    if (any(is.na(prior.beta$var))) stop("NA in prior.beta$var")
    if (any(prior.beta$var <= 0)) stop(paste("prior.beta$var must be higher than ", 0, sep=""))    
    Cbetavar <- as.numeric(prior.beta$var)
    names(Cbetavar) <- names(prior.beta$var) <- paste("beta", 1:lbeta, ".var", sep="")

    ##### put all together
    ##### -----------------------------------------------
    CpriorDouble_beta <- c(Cbetamean, Cbetavar)        
  }else{
    prior.beta <- list(mean=0, var=0)
    CpriorDouble_beta <- rep(0, 2)
  }  

  RET <- list(prior.beta        = prior.beta,
              CpriorDouble_beta = CpriorDouble_beta)
  return(RET)  
}
