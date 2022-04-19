\name{GLMM_MCMC}
\alias{GLMM_MCMC}
\alias{print.GLMM_MCMC}
\alias{print.GLMM_MCMClist}
\alias{C_GLMM_MCMC}
\alias{C_GLMM_PED}
\title{
  MCMC estimation of a (multivariate) generalized linear mixed model
  with a normal mixture in the distribution of random effects
}
\description{
  This function runs MCMC for a generalized linear mixed model
  with possibly several response variables and possibly normal
  mixtures in the distributions of random effects.
}  
\details{
  See accompanying papers (Komárek et al., 2010, Komárek and Komárková, 2013).  
}  
\usage{
GLMM_MCMC(y, dist = "gaussian", id, x, z, random.intercept,
     prior.alpha, init.alpha, init2.alpha,                      
     scale.b,     prior.b,    init.b,      init2.b,
     prior.eps,   init.eps,   init2.eps,
     nMCMC = c(burn = 10, keep = 10, thin = 1, info = 10),
     tuneMCMC = list(alpha = 1, b = 1),
     store = c(b = FALSE), PED = TRUE, keep.chains = TRUE,
     dens.zero = 1e-300, parallel = FALSE, cltype, silent = FALSE)

\method{print}{GLMM_MCMC}(x, \dots)

\method{print}{GLMM_MCMClist}(x, \dots)
}
\arguments{
  \item{y}{vector, matrix or data frame with responses. If \code{y} is vector then
     there is only one response in the model. If \code{y} is matrix or data frame then
     each column gives values of one response. Missing values are allowed.

     If there are several responses specified then 
     continuous responses must be put in the first columns and
     discrete responses in the subsequent columns.
  }
  \item{dist}{character (vector) which determines distribution (and a link function)
     for each response variable. Possible values are:
     \dQuote{gaussian} for gaussian (normal) distribution (with identity link),
     \dQuote{binomial(logit)} for binomial (0/1) distribution with a logit link.
     \dQuote{poisson(log)} for Poisson distribution with a log link.     
     Single value is recycled if necessary.
  }
  \item{id}{vector which determines longitudinally or otherwise dependent observations. If not given then it is assumed
    that there are no clusters and all observations of one response are
    independent.
  }
  \item{x}{matrix or a list of matrices with covariates (intercept not included) for fixed effects.
    If there is more than one response, this must always be a list. Note that intercept in included
    in all models. Use a character value \dQuote{empty} as a component of the list \code{x}
    if there are no covariates for a particular response.
  }
  \item{z}{matrix or a list of matrices with covariates (intercept not included) for random effects.
    If there is more than one response, this must always be a list. Note that random intercept
    is specified using the argument \code{random.intercept}.

    REMARK: For a particular response, matrices \code{x} and \code{z}
    may not have the same columns. That is, matrix \code{x} includes
    covariates which are not involved among random effects and matrix
    \code{z} includes covariates which are involved among random effects
    (and implicitely among fixed effects as well).
  }
  \item{random.intercept}{logical (vector) which determines for which responses random intercept
    should be included.
  }
  \item{prior.alpha}{a list which specifies prior distribution for fixed
    effects (not the means of random effects). The prior distribution is
    normal and the user can specify the mean and variances.
    The list \code{prior.alpha} can have the components listed below.
    \describe{
      \item{mean}{a vector with prior means, defaults to zeros.}
      \item{var}{a vector with prior variances, defaults to 10000 for
	all components.}      
    }      
  }
  \item{init.alpha}{a numeric vector with initial values of fixed effects
    (not the means of random effects) for the first chain. A sensible value is determined using the
    maximum-likelihood fits (using \code{\link[lme4]{lmer}} functions)
    and does not have to be given by the user.
  }
  \item{init2.alpha}{a numeric vector with initial values of fixed
    effects for the second chain.
  }  
  \item{scale.b}{a list specifying how to scale the random effects during
    the MCMC. A sensible value is determined using the
    maximum-likelihood fits (using \code{\link[lme4]{lmer}} functions)
    and does not have to be given by the user.

    If the user wishes to influence the shift and scale constants, these
    are given as components of the list \code{scale.b}. The components
    are named:
    \describe{
      \item{shift}{see vignette \href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{scale}{see vignette \href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
    }          
  }  
  \item{prior.b}{a list which specifies prior distribution for (shifted
    and scaled) random effects. The prior is in principle a normal
    mixture (being a simple normal distribution if we restrict the
    number of mixture components to be equal to one).
    
    The list \code{prior.b} can have the components listed below. Their
    meaning is analogous to the components of the same name of the
    argument \code{prior} of function \code{\link{NMixMCMC}} (see
    therein for details).
    \describe{
      \item{distribution}{a character string which specifies the assumed
        prior distribution for random effects. It can be either
	\dQuote{normal} (multivaruate normal - default) or \dQuote{MVT}
	(multivariate Student t distribution).}
      \item{priorK}{a character string which specifies the type of the
	prior for \eqn{K} (the number of mixture components).
      }
      \item{priormuQ}{a character string which specifies the type of the
        prior for mixture means and mixture variances.
      }
      \item{Kmax}{maximal number of mixture components.}
      \item{lambda}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{delta}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{xi}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{ce}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{D}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{zeta}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{gD}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{hD}{see vignette
	\href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
      \item{gdf}{shape parameter of the prior distribution for the
	degrees of freedom if the random effects are assumed to follow
	the MVT distribution}
      \item{hdf}{rate parameter of the prior distribution for the
	degrees of freedom if the random effects are assumed to follow
	the MVT distribution}
    }  
  }
  \item{init.b}{a list with initial values of the first chain for parameters related to the
    distribution of random effects and random effects
    themselves. Sensible initial values are determined by the function
    itself and do not have to be given by the user.
    \describe{
      \item{b}{}
      \item{K}{}
      \item{w}{}
      \item{mu}{}
      \item{Sigma}{}
      \item{Li}{}
      \item{gammaInv}{}
      \item{df}{}
      \item{r}{}                  
    }  
  }
  \item{init2.b}{a list with initial values of the second chain for parameters related to the
    distribution of random effects and random effects
    themselves.}
  \item{prior.eps}{a list specifying prior distributions for
     error terms for continuous responses. The list \code{prior.eps} can
     have the components listed below. For all components, a sensible
     value leading to weakly informative prior distribution can be
     determined by the function.
     \describe{
       \item{zeta}{see vignette
	 \href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
       \item{g}{see vignette
	 \href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}
       \item{h}{see vignette
	 \href{http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf}{PBCseq.pdf} for details}       
     }  
  }
  \item{init.eps}{a list with initial values of the first chain for parameters related to the
    distribution of error terms of continuous responses. The list
    \code{init.eps} can have the components listed below. For all
    components, a sensible value can be determined by the function.
    \describe{
      \item{sigma}{a numeric vector with the initial values for residual
        standard deviations for each continuous response.}
      \item{gammaInv}{a numeric vector with the initial values for the
	inverted components of the hyperparameter gamma for each
	continuous response.}
    }  
  }
  \item{init2.eps}{a list with initial values of the second chain for parameters related to the
    distribution of error terms of continuous responses.}  
  \item{nMCMC}{numeric vector of length 4 giving parameters of the MCMC
    simulation. Its components may be named (ordering is then unimportant) as:
    \describe{
      \item{burn}{length of the burn-in (after discarding the thinned
	values), can be equal to zero as well.}
      \item{keep}{length of the kept chains (after discarding the
	thinned values), must be positive.}
      \item{thin}{thinning interval, must be positive.}
      \item{info}{interval in which the progress information is printed
	on the screen.}      
    }
    In total \eqn{(M_{burn} + M_{keep}) \cdot M_{thin}}{(M[burn] +
      M[keep]) * M[thin]} MCMC scans are performed.
  }
  \item{tuneMCMC}{a list with tuning scale parameters for proposal
    distribution of fixed and random effects. It is used only when there
    are some discrete response profiles. The components of the list have
    the following meaning:
    \describe{
      \item{alpha}{scale parameters by which we multiply the proposal
	covariance matrix when updating the fixed effects pertaining to
	the discrete response profiles. There is one scale parameter for
        each DISCRETE profile. A single value is recycled if necessary.}
      \item{b}{a scale parameter by which we multiply the proposal
	covariance matrix when updating the random effects. It is used
	only when there are some discrete response profiles in the model.}
    }           
  }
  \item{store}{logical vector indicating whether the chains of
    parameters should be stored. Its components may be named (ordering
    is then unimportant) as:
    \describe{
      \item{b}{if \code{TRUE} then the sampled values of random effects
	are stored. Defaults to \code{FALSE}.}
    }  
  }
  \item{PED}{a logical value which indicates whether the penalized
    expected deviance (see Plummer, 2008 for more details)
    is to be computed (which requires two parallel
    chains).
  }  
  \item{keep.chains}{logical. If \code{FALSE}, only summary statistics
    are returned in the resulting object. This might be useful in the
    model searching step to save some memory.
  }
  \item{dens.zero}{a small value used instead of zero when computing
    deviance related quantities.
  }  
  \item{parallel}{a logical value which indicates whether parallel
    computation (based on a package \code{parallel})
    should be used when running two chains for the purpose of \code{PED}
    calculation.}
  \item{cltype}{optional argument applicable if \code{parallel} is
    \code{TRUE}. If \code{cltype} is given, it is passed as the
    \code{type} argument into the call to \code{\link[parallel]{makeCluster}}.}
  \item{silent}{a logical value indicating whether the information on
    the MCMC progress is to be supressed.}
  \item{\dots}{additional arguments passed to the default \code{print} method.}  
}
\value{
  An object of class \code{GLMM_MCMClist} (if \code{PED} argument is
  \code{TRUE}) or \code{GLMM_MCMC} (if \code{PED} argument is
  \code{FALSE}).
}
\section{Object of class GLMM_MCMC}{
  Object of class \code{GLMM_MCMC} can have the following
  components (some of them may be missing according to the context
  of the model):
  \describe{  %%% describe
  \item{iter}{index of the last iteration performed.}
  \item{nMCMC}{used value of the argument \code{nMCMC}.}
  \item{dist}{a character vector of length R corresponding to the
    \code{dist} argument.}  
  \item{R}{a two component vector giving the number of continuous
    responses and the number of discrete responses.}
  \item{p}{a numeric vector of length R giving the number of
    non-intercept alpha parameters for each response.}
  \item{q}{a numeric vector of length R giving the number of
    non-intercept random effects for each response.}
  \item{fixed.intercept}{a logical vector of length R which indicates
    inclusion of fixed intercept for each response.}
  \item{random.intercept}{a logical vector of length R which indicates
    inclusion of random intercept for each response.}  
  \item{lalpha}{length of the vector of fixed effects.}
  \item{dimb}{dimension of the distribution of random effects.}
  \item{prior.alpha}{a list containing the used value of the
    argument \code{prior.alpha}.}
  \item{prior.b}{a list containing the used value of the
    argument \code{prior.b}.}
  \item{prior.eps}{a list containing the used value of the
    argument \code{prior.eps}.}
  \item{init.alpha}{a numeric vector with the used value of the
    argument \code{init.alpha}.}
  \item{init.b}{a list containing the used value of the
    argument \code{init.b}.}
  \item{init.eps}{a list containing the used value of the
    argument \code{init.eps}.}
  \item{state.first.alpha}{a numeric vector with the first stored (after burn-in) value
    of fixed effects \eqn{\alpha}{alpha}.}
  \item{state.last.alpha}{a numeric vector with the last sampled value
    of fixed effects \eqn{\alpha}{alpha}. It can be used as argument
    \code{init.alpha} to restart MCMC.}
  \item{state.first.b}{a list with the first stored (after burn-in) values of parameters
    related to the distribution of random effects. It has components
    named \code{b}, \code{K}, \code{w}, \code{mu}, \code{Sigma}, \code{Li}, \code{Q},
    \code{gammaInv}, \code{r}.}
  \item{state.last.b}{a list with the last sampled values of parameters
    related to the distribution of random effects. It has components
    named \code{b}, \code{K}, \code{w}, \code{mu}, \code{Sigma}, \code{Li}, \code{Q},
    \code{gammaInv}, \code{r}. It can be used as argument
    \code{init.b} to restart MCMC.}
  \item{state.first.eps}{a list with the first stored (after burn-in) values of parameters
    related to the distribution of residuals of continuous responses. It
    has components named \code{sigma}, \code{gammaInv}.}
  \item{state.last.eps}{a list with the last sampled values of parameters
    related to the distribution of residuals of continuous responses. It
    has components named \code{sigma}, \code{gammaInv}. It can be used as argument
    \code{init.eps} to restart MCMC.}
  \item{prop.accept.alpha}{acceptance proportion from the
    Metropolis-Hastings algorithm for fixed effects
    (separately for each response type). Note that the acceptance
    proportion is equal to one for continuous responses since the Gibbs
    algorithm is used there.}
  \item{prop.accept.b}{acceptance proportion from the
    Metropolis-Hastings algorithm for random effects
    (separately for each cluster). Note that the acceptance
    proportion is equal to one for models with continuous responses only
    since the Gibbs algorithm is used there.}  
  \item{scale.b}{a list containing the used value of the argument
    \code{scale.b}.}
  \item{summ.Deviance}{a \code{data.frame} with posterior summary
    statistics for the deviance (approximated using the Laplacian
    approximation) and conditional (given random effects) devience.}
  \item{summ.alpha}{a \code{data.frame} with posterior summary statistics for fixed effects.}
  \item{summ.b.Mean}{a matrix with posterior summary statistics for
    means of random effects.}
  \item{summ.b.SDCorr}{a matrix with posterior summary statistics for
    standard deviations of random effects and correlations of each pair
    of random effects.}
  \item{summ.sigma_eps}{a matrix with posterior summary statistics for
    standard deviations of the error terms in the (mixed) models of
    continuous responses.}
  \item{poster.comp.prob_u}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob_u} is a matrix with \eqn{K} columns and \eqn{I}
    rows (\eqn{I} is the number of subjects defining the longitudinal
    profiles or correlated observations) with estimated posterior component probabilities
    -- posterior means of the components of the underlying 0/1
    allocation vector.

    \bold{WARNING:} By default, the labels of components are based on artificial
    identifiability constraints based on ordering of the mixture means
    in the first margin. Very often, such identifiability constraint is
    not satisfactory!
    
    %These can be used for possible clustering of the subjects based on
    %the longitudinal profiles.
  }
  \item{poster.comp.prob_b}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob_b} is a matrix with \eqn{K} columns and \eqn{I}
    rows (\eqn{I} is the number of subjects defining the longitudinal
    profiles or correlated observations)
    with estimated posterior component probabilities
    -- posterior mean over model parameters including random effects.

    \bold{WARNING:} By default, the labels of components are based on artificial
    identifiability constraints based on ordering of the mixture means
    in the first margin. Very often, such identifiability constraint is
    not satisfactory!
    
    %These can be used for possible clustering of the subjects based on
    %the longitudinal profiles.
  }    
  \item{freqK_b}{frequency table for the MCMC sample of the number of
    mixture components in the distribution of the random effects.}
  \item{propK_b}{posterior probabilities for the numbers of mixture
    components in the distribution of random effects.}
  \item{poster.mean.y}{a list with \code{data.frame}s, one
    \code{data.frame} per response profile. Each \code{data.frame}
    with columns labeled \code{id}, \code{observed},
    \code{fitted}, \code{stres},
    \code{eta.fixed} and \code{eta.random} holding
    identifier for clusters of grouped observations,
    observed values and
    posterior means for fitted values (response expectation given fixed and random effects),
    standardized residuals (derived from fitted values),
    fixed effect part of the linear predictor and the random effect part of
    the linear predictor. In each column, there are first all values for
    the first response, then all values for the second response etc.}
  \item{poster.mean.profile}{a \code{data.frame} with columns labeled
    \code{b1}, \ldots, \code{bq}, \code{Logpb}, \code{Cond.Deviance}, \code{Deviance} with
    posterior means of random effects for each cluster, posterior
    means of \eqn{\log\bigl\{p(\boldsymbol{b})\bigr\}}{log(p(b))},
    conditional deviances, i.e., minus twice the conditional (given
    random effects) log-likelihood for each cluster
    and GLMM deviances, i.e., minus twice the marginal (random effects
    integrated out) log-likelihoods for each cluster. The value of the
    marginal (random effects integrated out) log-likelihood at each MCMC
    iteration is obtained using the Laplacian approximation.}
  \item{poster.mean.w_b}{a numeric vector with posterior means of mixture
    weights after re-labeling. It is computed only if \eqn{K_b}{K[b]} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.mu_b}{a matrix with posterior means of mixture
    means after re-labeling. It is computed only if \eqn{K_b}{K[b]} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Q_b}{a list with posterior means of mixture inverse
    variances after re-labeling. It is computed only if \eqn{K_b}{K[b]} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Sigma_b}{a list with posterior means of mixture
    variances after re-labeling. It is computed only if \eqn{K_b}{K[b]} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Li_b}{a list with posterior means of Cholesky
    decompositions of mixture inverse
    variances after re-labeling. It is computed only if \eqn{K_b}{K[b]} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{Deviance}{numeric vector with a chain for the GLMM deviances,
    i.e., twice the marginal (random effects integrated out)
    log-likelihoods of the GLMM. The marginal log-likelihood is obtained
    using the Laplacian approximation at each iteration of MCMC.
  }
  \item{Cond.Deviance}{numeric vector with a chain for the conditional deviances,
    i.e., twice the conditional (given random effects)
    log-likelihoods.
  }  
  \item{K_b}{numeric vector with a chain for \eqn{K_b}{K[b]} (number of
    mixture components in the distribution of random effects).}
  \item{w_b}{numeric vector or matrix with a chain for \eqn{w_b}{w[b]} (mixture
    weights for the distribution of random effects). It is a matrix with
    \eqn{K_b}{K[b]} columns when \eqn{K_b}{K[b]} is fixed.
    Otherwise, it is a vector with weights put sequentially after each other.}
  \item{mu_b}{numeric vector or matrix with a chain for \eqn{\mu_b}{mu[b]} (mixture
    means for the distribution of random effects). It is a matrix with
    \eqn{dimb\cdot K_b}{dimb*K[b]} columns when \eqn{K_b}{K[b]} is
    fixed. Otherwise, it is a vector with means put sequentially after
    each other.}
  \item{Q_b}{numeric vector or matrix with a chain for lower triangles of \eqn{\boldsymbol{Q}_b}{Q[b]} (mixture
    inverse variances for the distribution of random effects). It is a matrix with \eqn{\frac{dimb(dimb+1)}{2}\cdot K_b}{(dimb*(dimb+1)2)*K[b]}
    columns when \eqn{K_b}{K[b]} is fixed. Otherwise, it is a vector with lower
    triangles of \eqn{\boldsymbol{Q}_b}{Q[b]} matrices put sequentially after each other.}
  \item{Sigma_b}{numeric vector or matrix with a chain for lower triangles of \eqn{\Sigma_b}{Sigma[b]} (mixture
    variances for the distribution of random effects). It is a matrix with \eqn{\frac{dimb(dimb+1)}{2}\cdot K_b}{(dimb*(dimb+1)2)*K[b]}
    columns when \eqn{K_b}{K[b]} is fixed. Otherwise, it is a vector with lower
    triangles of \eqn{\Sigma_b}{Sigma[b]} matrices put sequentially after each other.}
  \item{Li_b}{numeric vector or matrix with a chain for lower triangles of
    Cholesky decompositions of \eqn{\boldsymbol{Q}_b}{Q[b]} matrices.
    It is a matrix with \eqn{\frac{dimb(dimb+1)}{2}\cdot K_b}{(dimb*(dimb+1)2)*K[b]}
    columns when \eqn{K_b}{K[b]} is fixed. Otherwise, it is a vector with lower
    triangles put sequentially after each other.}  
  \item{gammaInv_b}{matrix with \eqn{dimb} columns with a chain for inverses
    of the hyperparameter \eqn{\boldsymbol{\gamma}_b}{gamma[b]}.}
  \item{order_b}{numeric vector or matrix with order indeces of mixture
    components in the distribution of random effects related to artificial identifiability constraint
    defined by ordering of the first component of the mixture means.

    It is a matrix with \eqn{K_b}{K[b]} columns when \eqn{K_b}{K[b]} is
    fixed. Otherwise it is a vector with orders put sequentially after
    each other.}
  \item{rank_b}{numeric vector or matrix with rank indeces of mixture
    components in the distribution of random effects related to artificial identifiability constraint
    defined by ordering of the first component of the mixture means.

    It is a matrix with \eqn{K_b}{K[b]} columns when \eqn{K_b}{K[b]} is
    fixed. Otherwise it is a vector with ranks put sequentially after
    each other.}
  \item{mixture_b}{\code{data.frame} with columns labeled
    \code{b.Mean.*}, \code{b.SD.*}, \code{b.Corr.*.*}
    containing the chains for the means, standard deviations and correlations of the
    distribution of the random effects based on a normal mixture at each
    iteration.}
  \item{b}{a matrix with the MCMC chains for random effects. It is
    included only if \code{store[b]} is \code{TRUE}.}
  \item{alpha}{numeric vector or matrix with the MCMC chain(s) for fixed effects.}
  \item{sigma_eps}{numeric vector or matrix with the MCMC chain(s) for
    standard deviations of the error terms in the (mixed) models for
    continuous responses.}
  \item{gammaInv_eps}{matrix with \eqn{dimb} columns with MCMC chain(s) for inverses
    of the hyperparameter \eqn{\boldsymbol{\gamma}_b}{gamma[b]}.}
  \item{relabel_b}{a list which specifies the algorithm used to re-label
    the MCMC output to compute \code{order_b}, \code{rank_b}, \code{poster.comp.prob_u},
    \code{poster.comp.prob_b}, \code{poster.mean.w_b},
    \code{poster.mean.mu_b}, \code{poster.mean.Q_b},
    \code{poster.mean.Sigma_b}, \code{poster.mean.Li_b}. 
  }
  \item{Cpar}{a list with components useful to call underlying C++ functions (not
    interesting for ordinary users).}
  }
}
\section{Object of class GLMM_MCMClist}{
  Object of class \code{NMixMCMClist} is the list having two components
  of class \code{NMixMCMC} representing two parallel chains and
  additionally the following components:
  \describe{
  \item{PED}{values of penalized expected deviance and related
    quantities. It is a vector with five components: \code{D.expect} \eqn{=}
    estimated expected deviance, where the estimate is based on two
    parallel chains; \code{popt} \eqn{=} estimated penalty, where the
    estimate is based on simple MCMC average based on two parallel
    chains; \code{PED} \eqn{=} estimated penalized expected deviance
    \eqn{=} \code{D.expect} \eqn{+} \code{popt}; \code{wpopt} \eqn{=}
    estimated penalty, where the estimate is based on weighted MCMC average
    (through importance sampling) based on two parallel chains;
    \code{wPED} \eqn{=} estimated penalized expected deviance \eqn{=}
    \code{D.expect} \eqn{+} \code{wpopt}.}
  \item{D}{posterior mean of the deviance for each subject.}
  \item{popt}{contributions to the unweighted penalty from each subject.}
  \item{wpopt}{contributions to the weighted penalty from each subject.}
  \item{inv.D}{for each subject, number of iterations (in both chains), where the
    deviance was in fact equal to infinity (when the corresponding
    density was lower than \code{dens.zero}) and was not taken into account when
    computing \code{D.expect}.}
  \item{inv.popt}{for each subject, number of iterations, where the
    penalty was in fact equal to infinity and was not taken into account
    when computing \code{popt}.
  }
  \item{inv.wpopt}{for each subject, number of iterations, where the
    importance sampling weight was in fact equal to infinity and was not taken into account
    when computing \code{wpopt}.}
  \item{sumISw}{for each subject, sum of importance sampling
    weights.}
  \item{Deviance1}{sampled value of the observed data deviance from
    chain 1}
  \item{Deviance2}{sampled values of the obserbed data deviance from
    chain 2}
  \item{Deviance_repl1_ch1}{sampled values of the deviance of data
    replicated according to the chain 1 evaluated under the parameters
    from chain 1}
  \item{Deviance_repl1_ch2}{sampled values of the deviance of data
    replicated according to the chain 1 evaluated under the parameters
    from chain 2}
  \item{Deviance_repl2_ch1}{sampled values of the deviance of data
    replicated according to the chain 2 evaluated under the parameters
    from chain 1}
  \item{Deviance_repl2_ch2}{sampled values of the deviance of data
    replicated according to the chain 2 evaluated under the parameters
    from chain 2}  
  }
}  
\references{
  Komárek, A. and Komárková, L. (2013).
  Clustering for multivariate continuous and discrete longitudinal data.
  \emph{The Annals of Applied Statistics}, \bold{7}(1), 177--200.

  Komárek, A. and Komárková, L. (2014).
  Capabilities of R package mixAK for clustering based on multivariate continuous and discrete longitudinal data.
  \emph{Journal of Statistical Software}, \bold{59}(12), 1--38.
  \doi{10.18637/jss.v059.i12}.
  
  Komárek, A., Hansen, B. E., Kuiper, E. M. M., van Buuren, H. R., and Lesaffre, E. (2010).
  Discriminant analysis using a multivariate linear mixed model with a normal mixture in the random effects distribution.
  \emph{Statistics in Medicine}, \bold{29}(30), 3267--3283.

  Plummer, M. (2008).
  Penalized loss functions for Bayesian model comparison.
  \emph{Biostatistics}, \bold{9}(3), 523--539.  
}  
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
## See also additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## - files http://www.karlin.mff.cuni.cz/~komarek/software/mixAK/PBCseq.pdf,
##         PBCseq.R
## ==============================================
}
\keyword{models}
\keyword{multivariate}

