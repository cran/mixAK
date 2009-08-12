\name{GLMM_MCMC}
\alias{GLMM_MCMC}
\alias{print.GLMM_MCMC}
\title{
  MCMC estimation of generalized linear mixed model
  with mixtures in the distributions.
}
\description{
  THIS FUNCTION IS BEING DEVELOPED AND ORDINARY USERS
  ARE NOT RECOMMENDED TO USE IT.
  
  This function runs MCMC for a generalized linear mixed model
  with possibly several response variables and possibly normal
  mixtures in the distributions of random effects.
}  
\details{
  I HOPE TO WRITE A PAPER DESCRIBING ALL THE DETAILS.
}  
\usage{
GLMM_MCMC(y, dist="gaussian", id, x, z, random.intercept,
          prior.beta, init.beta,                      
          scale.b,    prior.b,   init.b,
          prior.eps,  init.eps,
          nMCMC=c(burn=10, keep=10, thin=1, info=10),
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
     \dQuote{gaussian} for gaussian (normal) distribution (with identity
     link),
     \dQuote{binomial(logit)} for binomial (0/1) distribution with a logit link.
     Single value is recycled if necessary.
  }
  \item{id}{vector which determines clustered observations. If not given then it is assumed
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
  \item{R}{a two component vector giving the number of continuous
    responses and the number of discrete responses.}
  \item{dist}{a character vector of length R corresponding to the
    \code{dist} argument.}
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
  \item{state.beta}{a numeric vector with the last sampled value
    of fixed effects \eqn{\beta}{beta}. It can be used as argument
    \code{init.beta} to restart MCMC.}
  \item{state.b}{a list with the last sampled values of parameters
    related to the distribution of random effects. It has components
    named \code{b}, \code{K}, \code{w}, \code{mu}, \code{Sigma}, \code{Li}, \code{Q},
    \code{gammaInv}, \code{r}. It can be used as argument
    \code{init.b} to restart MCMC.}
  \item{state.eps}{a list with the last sampled values of parameters
    related to the distribution of residuals of continuous responses. It
    has components named \code{sigma}, \code{gammaInv}. It can be used as argument
    \code{init.eps} to restart MCMC.}
  \item{scale.b}{a list containing the used value of the argument
    \code{scale.b}.}
  \item{poster.mean.eta}{a \code{data.frame} with columns labeled
    \code{fixed} and \code{random} holding posterior means for fixed
    effect part of the linear predictor and the random effect part of
    the linear predictor. In each column, there are first all values for
    the first response, then all values for the second response etc.}
  \item{poster.mean.cluster}{a \code{data.frame} with columns labeled
    \code{b1}, \ldots, \code{bq}, \code{LogL}, \code{Logpb} with
    posterior means of random effects for each cluster and posterior
    means of \eqn{\log(L)}{log(L)} (log-likelihood given random effects)
    and \eqn{\log\bigl\{p(\boldsymbol{b})\bigr\}}{log(p(b))} for each cluster.}
  \item{poster.mean.w_b}{ADD DESCRIPTION}
  \item{poster.mean.mu_b}{ADD DESCRIPTION}
  \item{poster.mean.Q_b}{ADD DESCRIPTION}
  \item{poster.mean.Sigma_b}{ADD DESCRIPTION}
  \item{poster.mean.Li_b}{ADD DESCRIPTION}
  \item{summ.beta}{ADD DESCRIPTION}
  \item{summ.b.Mean}{ADD DESCRIPTION}
  \item{summ.b.SDCorr}{ADD DESCRIPTION}
  \item{summ.sigma_eps}{ADD DESCRIPTION}        
  \item{freqK_b}{ADD DESCRIPTION}
  \item{propK_b}{ADD DESCRIPTION}  
  \item{K_b}{ADD DESCRIPTION}
  \item{w_b}{ADD DESCRIPTION}
  \item{mu_b}{ADD DESCRIPTION}
  \item{Li_b}{ADD DESCRIPTION}
  \item{Q_b}{ADD DESCRIPTION}
  \item{Sigma_b}{ADD DESCRIPTION}
  \item{gammaInv_b}{ADD DESCRIPTION}
  \item{order_b}{ADD DESCRIPTION}
  \item{rank_b}{ADD DESCRIPTION}
  \item{mixture_b}{ADD DESCRIPTION}
  \item{b}{ADD DESCRIPTION}
  \item{beta}{ADD DESCRIPTION}
  \item{sigma_eps}{ADD DESCRIPTION}
  \item{gammaInv_eps}{ADD DESCRIPTION}
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.} (20XX).
  PAPER TO BE WRITTEN.
  \emph{JOURNAL}, \bold{XX}, XXX-XXX.
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

