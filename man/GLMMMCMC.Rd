\name{GLMM_MCMC}
\alias{GLMM_MCMC}
\alias{print.GLMM_MCMC}
\title{
  MCMC estimation of generalized linear mixed model
  with mixtures in the distributions.
}
\description{
  This function runs MCMC for a generalized linear mixed model
  with possibly several response variables and possibly normal
  mixtures in the distributions of random effects.
}  
\details{
  See accompanying paper (\eqn{\mbox{Kom\'arek}}{Kom&#225rek} et al., 2010).  
}  
\usage{
GLMM_MCMC(y, dist="gaussian", id, x, z, random.intercept,
          prior.beta, init.beta,                      
          scale.b,    prior.b,   init.b,
          prior.eps,  init.eps,
          nMCMC=c(burn=10, keep=10, thin=1, info=10),
          tuneMCMC=list(beta=1, b=1),
          store=c(b=FALSE), keep.chains=TRUE)

\method{print}{GLMM_MCMC}(x, \dots)
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
  \item{prior.beta}{a list which specifies prior distribution for fixed
    effects (not the means of random effects). The prior distribution is
    normal and the user can specify the mean and variances.
    The list \code{prior.b} can have the components listed below.
    \describe{
      \item{mean}{a vector with prior means, defaults to zeros.}
      \item{var}{a vector with prior variances, defaults to 10000 for
	all components.}      
    }      
  }
  \item{init.beta}{a numeric vector with initial values of fixed effects
    (not the means of random effects). A sensible value is determined using the
    maximum-likelihood fits (using \code{\link[lme4]{lmer}} functions)
    and does not have to be given by the user.
  }
  \item{scale.b}{a list specifying how to scale the random effects during
    the MCMC. A sensible value is determined using the
    maximum-likelihood fits (using \code{\link[lme4]{lmer}} functions)
    and does not have to be given by the user.

    If the user wishes to influence the shift and scale constants, these
    are given as components of the list \code{scale.b}. The components
    are named:
    \describe{
      \item{shift}{ADD DESCRIPTION}
      \item{scale}{ADD DESCRIPTION}
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
      \item{priorK}{a character string which specifies the type of the
	prior for \eqn{K} (the number of mixture components).
      }
      \item{priormuQ}{a character string which specifies the type of the
        prior for mixture means and mixture variances.
      }
      \item{Kmax}{maximal number of mixture components.}
      \item{lambda}{ADD DESCRIPTION}
      \item{delta}{ADD DESCRIPTION}
      \item{xi}{ADD DESCRIPTION}
      \item{ce}{ADD DESCRIPTION}
      \item{D}{ADD DESCRIPTION}
      \item{zeta}{ADD DESCRIPTION}
      \item{g}{ADD DESCRIPTION}
      \item{h}{ADD DESCRIPTION}            
    }  
  }
  \item{init.b}{a list with initial values for parameters related to the
    distribution of random effects and random effects
    themselves. Sensible initial values are determined by the function
    itself and do not have to be given by the user.
  }
  \item{prior.eps}{a list specifying prior distributions for
     error terms for continuous responses. The list \code{prior.eps} can
     have the components listed below. For all components, a sensible
     value leading to weakly informative prior distribution can be
     determined by the function.
     \describe{
       \item{zeta}{ADD DESCRIPTION}
       \item{g}{ADD DESCRIPTION}
       \item{h}{ADD DESCRIPTION}       
     }  
  }
  \item{init.eps}{a list with initial values for parameters related to the
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
      \item{beta}{scale parameters by which we multiply the proposal
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
  \item{keep.chains}{logical. If \code{FALSE}, only summary statistics
    are returned in the resulting object. This might be useful in the
    model searching step to save some memory.
  }
  \item{\dots}{additional arguments passed to the default \code{print} method.}  
}
\value{
  An object of class \code{GLMM_MCMC}. It can have the following
  components (some of them may be missing according to the context
  of the model):
  \item{iter}{index of the last iteration performed.}
  \item{nMCMC}{used value of the argument \code{nMCMC}.}
  \item{dist}{a character vector of length R corresponding to the
    \code{dist} argument.}  
  \item{R}{a two component vector giving the number of continuous
    responses and the number of discrete responses.}
  \item{p}{a numeric vector of length R giving the number of
    non-intercept beta parameters for each response.}
  \item{q}{a numeric vector of length R giving the number of
    non-intercept random effects for each response.}
  \item{fixed.intercept}{a logical vector of length R which indicates
    inclusion of fixed intercept for each response.}
  \item{random.intercept}{a logical vector of length R which indicates
    inclusion of random intercept for each response.}  
  \item{lbeta}{length of the vector of fixed effects.}
  \item{dimb}{dimension of the distribution of random effects.}
  \item{prior.beta}{a list containing the used value of the
    argument \code{prior.beta}.}
  \item{prior.b}{a list containing the used value of the
    argument \code{prior.b}.}
  \item{prior.eps}{a list containing the used value of the
    argument \code{prior.eps}.}
  \item{init.beta}{a numeric vector with the used value of the
    argument \code{init.beta}.}
  \item{init.b}{a list containing the used value of the
    argument \code{init.b}.}
  \item{init.eps}{a list containing the used value of the
    argument \code{init.eps}.}
  \item{state.first.beta}{a numeric vector with the first stored (after burn-in) value
    of fixed effects \eqn{\beta}{beta}.}
  \item{state.last.beta}{a numeric vector with the last sampled value
    of fixed effects \eqn{\beta}{beta}. It can be used as argument
    \code{init.beta} to restart MCMC.}
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
  \item{prop.accept.beta}{acceptance proportion from the
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
  \item{summ.beta}{a \code{data.frame} with posterior summary statistics for fixed effects.}
  \item{summ.b.Mean}{a matrix with posterior summary statistics for
    means of random effects.}
  \item{summ.b.SDCorr}{a matrix with posterior summary statistics for
    standard deviations of random effects and correlations of each pair
    of random effects.}
  \item{summ.sigma_eps}{a matrix with posterior summary statistics for
    standard deviations of the error terms in the (mixed) models of
    continuous responses.}
  \item{poster.comp.prob1}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob1} is a matrix with \eqn{K} columns and \eqn{I}
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
  \item{poster.comp.prob2}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob2} is a matrix with \eqn{K} columns and \eqn{I}
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
  \item{beta}{numeric vector or matrix with the MCMC chain(s) for fixed effects.}
  \item{sigma_eps}{numeric vector or matrix with the MCMC chain(s) for
    standard deviations of the error terms in the (mixed) models for
    continuous responses.}
  \item{gammaInv_eps}{matrix with \eqn{dimb} columns with MCMC chain(s) for inverses
    of the hyperparameter \eqn{\boldsymbol{\gamma}_b}{gamma[b]}.}
  \item{relabel_b}{a list which specifies the algorithm used to re-label
    the MCMC output to compute \code{order_b}, \code{rank_b}, \code{poster.comp.prob1},
    \code{poster.comp.prob2}, \code{poster.mean.w_b},
    \code{poster.mean.mu_b}, \code{poster.mean.Q_b},
    \code{poster.mean.Sigma_b}, \code{poster.mean.Li_b}. 
  }
  \item{Cpar}{a list with components useful to call underlying C++ functions (not
    interesting for ordinary users).
  }  
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.}, Hansen, B. E., Kuiper,
  E. M. M., van Buuren, H. R., and Lesaffre, E. (2010).
  Discriminant analysis using a multivariate linear mixed model with a
  normal mixture in the random effects distribution.
  \emph{Statistics in Medicine}. To appear. %\bold{XX}, XXX-XXX.
}  
\seealso{
  \code{\link{NMixMCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### WILL BE ADDED.
}
\keyword{models}
\keyword{multivariate}

