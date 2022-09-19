\name{NMixMCMC}
\alias{NMixMCMC}
\alias{print.NMixMCMC}
\alias{print.NMixMCMClist}
\alias{C_NMix_MCMC}
\alias{C_NMix_PED}
\title{
  MCMC estimation of (multivariate) normal mixtures with possibly
  censored data.
}
\description{
  This function runs MCMC for a model in which unknown density is
  specified as a normal mixture with either known or unknown number of
  components. With a prespecified number of components, MCMC is
  implemented through Gibbs sampling (see Diebolt and Robert, 1994) and
  dimension of the data can be arbitrary.
  With unknown number of components, currently only univariate case is
  implemented using the reversible jump MCMC (Richardson and Green,
  1997).
  %In a future, there are plans to implement the reversible jump
  %MCMC in a higher dimension using the algorithm presented by
  %Papageorgiou and Dellaportas (2006).
  
  Further, the data are allowed to be censored in which case additional
  Gibbs step is used within the MCMC algorithm 
}
\details{
  See accompanying paper (Komárek, 2009).
  In the rest of the helpfile,
  the same notation is used as in the paper, namely, \eqn{n} denotes the number of
  observations, \eqn{p} is dimension of the data, \eqn{K} is the number
  of mixture components,
  \eqn{w_1,\dots,w_K}{w[1],\ldots,w[K]} are mixture weights, 
  \eqn{\boldsymbol{\mu}_1,\dots,\boldsymbol{\mu}_K}{mu[1],\ldots,mu[K]}
  are mixture means,
  \eqn{\boldsymbol{\Sigma}_1,\dots,\boldsymbol{\Sigma}_K}{Sigma_1,\ldots,Sigma_K}
  are mixture variance-covariance matrices,
  \eqn{\boldsymbol{Q}_1,\dots,\boldsymbol{Q}_K}{Q_1,\ldots,Q_K} are
  their inverses.

  For the data
  \eqn{\boldsymbol{y}_1,\dots,\boldsymbol{y}_n}{y[1],\ldots,y[n]} the
  following \eqn{g_y(\boldsymbol{y})}{g(y)} density is assumed
  \deqn{
    g_y(\boldsymbol{y}) = |\boldsymbol{S}|^{-1} \sum_{j=1}^K w_j
    \varphi\bigl(\boldsymbol{S}^{-1}(\boldsymbol{y} - \boldsymbol{m}\,|\,\boldsymbol{\mu}_j,\,\boldsymbol{\Sigma}_j)\bigr),
  }{%
    g(y) = |S|^(-1) * sum[j=1]^K w[j] phi(S^(-1)*(y - m) | mu[j], Sigma[j]),
  }
  where
  \eqn{\varphi(\cdot\,|\,\boldsymbol{\mu},\,\boldsymbol{\Sigma})}{phi(. | mu, Sigma)} denotes a density
  of the (multivariate) normal distribution
  with mean \eqn{\boldsymbol{\mu}}{mu} and a~variance-covariance matrix \eqn{\boldsymbol{\Sigma}}{\Sigma}.
  Finally, \eqn{\boldsymbol{S}}{S} is a pre-specified diagonal scale matrix and
  \eqn{\boldsymbol{m}}{m} is a pre-specified shift vector. Sometimes, by
  setting \eqn{\boldsymbol{m}}{m} to sample means of components of
  \eqn{\boldsymbol{y}}{y} and diagonal of \eqn{\boldsymbol{S}}{S} to
  sample standard deviations of \eqn{\boldsymbol{y}}{y} (considerable)
  improvement of the MCMC algorithm is achieved.   
}
\usage{
NMixMCMC(y0, y1, censor, x_w, scale, prior,
         init, init2, RJMCMC,
         nMCMC = c(burn = 10, keep = 10, thin = 1, info = 10),
         PED, keep.chains = TRUE, onlyInit = FALSE, dens.zero = 1e-300,
         parallel = FALSE, cltype)

\method{print}{NMixMCMC}(x, dic, \dots)

\method{print}{NMixMCMClist}(x, ped, dic, \dots)
}
\arguments{
  \item{y0}{numeric vector of length \eqn{n} or \eqn{n\times p}{n x p}
    matrix with observed data. It contains exactly observed,
    right-censored, left-censored data and lower limits for
    interval-censored data.}
  \item{y1}{numeric vector of length \eqn{n} or \eqn{n\times p}{n x p}
    matrix with upper limits for interval-censored data. Elements
    corresponding to exactly observed, right-censored or left-censored
    data are ignored and can be filled arbitrarily (by
    \code{NA}'s) as well.

    It does not have to be supplied if there are no interval-censored data.
  }
  \item{censor}{numeric vector of length \eqn{n} or \eqn{n\times p}{n x p}
    matrix with censoring indicators. The following values indicate:

    \describe{
      \item{0}{right-censored observation,}
      \item{1}{exactly observed value,}
      \item{2}{left-censored observation,}
      \item{3}{interval-censored observation.}      
    }

    If it is not supplied then it is assumed that all values are exactly
    observed.
  }
  \item{x_w}{optional vector providing a categorical covariate that may
    influence the mixture weights. Internally, it is converted into a
    \code{factor}.

    Added in version 4.0 (03/2015).
  }
  \item{scale}{a list specifying how to scale the data before running
    MCMC. It should have two components:

    \describe{
      \item{shift}{a vector of length 1 or \eqn{p} specifying shift
	vector \eqn{\boldsymbol{m}}{m},}
      \item{scale}{a vector of length 1 or \eqn{p} specifying diagonal of the
        scaling matrix \eqn{\boldsymbol{S}}{S}.}
    }  

    If there is no censoring, and argument \code{scale} is missing
    then the data are scaled to have zero mean and unit variances, i.e.,
    \code{scale(y0)} is used for MCMC. In the case there is censoring
    and \code{scale} is missing,
    \code{scale$shift} is taken to be a sample mean of \code{init$y} and
    \code{scale$scale} are sample standard deviations of columns of \code{init$y}.
    
    If you do not wish to scale the data before running MCMC, specify
    \code{scale=list(shift=0, scale=1)}.
  }
  \item{prior}{a list with the parameters of the prior distribution. It
    should have the following components (for some of them,
    the program can assign default values and the user does not have to
    specify them if he/she wishes to use the defaults):
    \describe{
      \item{priorK}{a character string which specifies the type of the
	prior for \eqn{K} (the number of mixture components). It should
	have one of the following values:

        \eqn{\mbox{\hspace{6in}}}{\cr}
	\dQuote{\code{fixed}}\eqn{\mbox{\hspace{6in}}}{\cr}
	Number of mixture components is assumed to be fixed to
	\eqn{K_{max}}{K[max]}. This is a \bold{default} value.

        \eqn{\mbox{\hspace{6in}}}{\cr}
	\dQuote{\code{uniform}}\eqn{\mbox{\hspace{6in}}}{\cr}
        A priori \eqn{K \sim \mbox{Unif}\{1,\dots,K_{max}\}.}{K ~ Unif{1,\ldots,K[max]}.}
		
        \eqn{\mbox{\hspace{6in}}}{\cr}
	\dQuote{\code{tpoisson}}\eqn{\mbox{\hspace{6in}}}{\cr}
	A priori \eqn{K \sim
	  \mbox{truncated-Poiss}(\lambda,\,K_{max}).}{K ~
	  truncated-Poiss(lambda, K[max]).}
      }
      \item{priormuQ}{a character string which specifies the type of the
        prior for \eqn{\boldsymbol{\mu}_1,\dots,\boldsymbol{\mu}_{K_{max}}}{mu[1], \ldots,
	mu[K[max]]} (mixture means) and
        \eqn{\boldsymbol{Q}_1,\dots,\boldsymbol{Q}_{K_{max}}}{Q[1], \ldots,
        Q[K[max]]} (inverted mixture covariance matrices). It should
        have one of the following values:

        \eqn{\mbox{\hspace{6in}}}{\cr}
	\dQuote{\code{independentC}}\eqn{\mbox{\hspace{6in}}}{\cr}
        \eqn{\equiv}{=} independent conjugate prior (this is a
	\bold{default} value).
	That is, a priori
	\deqn{
          (\boldsymbol{\mu}_j,\, \boldsymbol{Q}_j) \sim
	  \mbox{N}(\boldsymbol{\xi}_j,\,\boldsymbol{D}_j)
	  \times
	  \mbox{Wishart}(\zeta,\,\boldsymbol{\Xi})
	}{
          (mu[j], Q[j]) ~ N(xi[j], D[j]) * Wishart(zeta, Xi) 
	}  
        independently for \eqn{j=1,\ldots,K}, where
	normal means
	\eqn{\boldsymbol{\xi}_1,\dots,\boldsymbol{\xi}_K}{xi[1],\ldots,xi[K]},
	normal variances
	\eqn{\boldsymbol{D}_1,\dots,\boldsymbol{D}_K}{D[1],\ldots,D[K]},
	and Wishart degrees of freedom \eqn{\zeta}{zeta} are specified
	further as \code{xi}, \code{D}, \code{zeta} components of the
	list \code{prior}.
	
        \eqn{\mbox{\hspace{6in}}}{\cr}
	\dQuote{\code{naturalC}}\eqn{\mbox{\hspace{6in}}}{\cr}
        \eqn{\equiv}{=} natural conjugate prior. That is, a priori
	\deqn{
          (\boldsymbol{\mu}_j,\, \boldsymbol{Q}_j) \sim
	  \mbox{N}(\boldsymbol{\xi}_j,\,c_j^{-1}\boldsymbol{Q}_j^{-1})
	  \times
	  \mbox{Wishart}(\zeta,\,\boldsymbol{\Xi})
	}{
          (mu[j], Q[j]) ~ N(xi[j], (c[j]Q[j])^(-1)) * Wishart(zeta, Xi) 
	}  
        independently for \eqn{j=1,\ldots,K}, where
	normal means
	\eqn{\boldsymbol{\xi}_1,\dots,\boldsymbol{\xi}_K}{xi[1],\ldots,xi[K]},
	precisions
	\eqn{c_1,\dots,c_K}{c[1],\ldots,c[K]},
	and Wishart degrees of freedom \eqn{\zeta}{zeta} are specified
	further as \code{xi}, \code{ce}, \code{zeta} components of the
	list \code{prior}.

        \eqn{\mbox{\hspace{6in}}}{\cr}
        For both, independent conjugate and natural conjugate prior, the
	Wishart scale matrix \eqn{\boldsymbol{\Xi}}{Xi} is assumed to be
	diagonal with
	\eqn{\gamma_1,\dots,\gamma_p}{gamma[1],\ldots,gamma[p]} on a
	diagonal. For \eqn{\gamma_j^{-1}}{gamma[j]^(-1)} \eqn{(j=1,\ldots,K)} additional
	gamma hyperprior \eqn{\mbox{G}(g_j,\,h_j)}{G(g[j], h[j])} is assumed.
	Values of \eqn{g_1,\dots,g_p}{g[1],\ldots,g[p]} and
	\eqn{h_1,\dots,h_p}{h[1],\ldots,h[p]} are further specified as
	\code{g} and \code{h} components of the \code{prior} list.
      }
      \item{Kmax}{maximal number of mixture components
	\eqn{K_{max}}{K[max]}. It must \bold{always be specified} by the user.}
      \item{lambda}{parameter \eqn{\lambda}{lambda} for the truncated
	Poisson prior on \eqn{K}. It must be positive and must \bold{always be specified} if
	\code{priorK} is \dQuote{\code{tpoisson}}.}
      \item{delta}{parameter \eqn{\delta}{delta} for the Dirichlet prior
        on the mixture weights \eqn{w_1,\dots,w_K.}{w[1],\ldots,w[K].}
        It must be positive. Its \bold{default} value is 1.}
      \item{xi}{a numeric value, vector or matrix which specifies
	\eqn{\boldsymbol{\xi}_1, \dots, \boldsymbol{\xi}_{K_{max}}}{xi[1],
	  \ldots, xi[K[max]]} (prior means for the mixture means
	  \eqn{\boldsymbol{\mu}_1,\dots,\boldsymbol{\mu}_{K_{max}}}{mu[1], \ldots,
	    mu[K[max]]}). \bold{Default} value is a matrix
	  \eqn{K_{max}\times p}{K[max] x p} with midpoints of columns
	  of \code{init$y} in rows which follows Richardson and Green (1997).

          \eqn{\mbox{\hspace{6in}}}{\cr}	  
          If \eqn{p=1} and \code{xi}\eqn{=\xi}{=xi} is a single value
	  then \eqn{\xi_1=\cdots=\xi_{K_{max}} =
	    \xi.}{xi[1]=\ldots=xi[K[max]] = xi.}
	  
          \eqn{\mbox{\hspace{6in}}}{\cr}
	  If \eqn{p=1} and \code{xi}\eqn{=\boldsymbol{\xi}} is a vector of length
	  \eqn{K_{max}}{K[max]} then the \eqn{j}-th element of \code{xi}
	  gives \eqn{\xi_j}{xi[j]} \eqn{(j=1,\dots,K_{max}).}{(j=1,\ldots,K[max]).}
	  
          \eqn{\mbox{\hspace{6in}}}{\cr}	  
          If \eqn{p>1} and \code{xi}\eqn{=\boldsymbol{\xi}}{=xi} is a vector of length \eqn{p}
	  then \eqn{\boldsymbol{\xi}_1=\cdots=\boldsymbol{\xi}_{K_{max}} =
	    \boldsymbol{\xi}.}{xi[1]=\ldots=xi[K[max]] = xi.}

          \eqn{\mbox{\hspace{6in}}}{\cr}
	  If \eqn{p>1} and \code{xi} is a \eqn{K_{max} \times p}{K[max]
	    x p} matrix then the \eqn{j}-th row of \code{xi}
	  gives \eqn{\boldsymbol{xi}_j}{xi[j]} \eqn{(j=1,\dots,K_{max}).}{(j=1,\ldots,K[max]).}	  
      } 
      \item{ce}{a numeric value or vector which specifies prior
	precision parameters
	\eqn{c_1,\dots,c_{K_{max}}}{c[1],\ldots,c[K[max]]} for the
	mixture means \eqn{\boldsymbol{\mu}_1,\dots,\boldsymbol{\mu}_{K_{max}}}{mu[1], \ldots,
          mu[K[max]]} when \code{priormuQ} is \dQuote{\code{naturalC}}.
        Its \bold{default} value is a vector of ones which follows
	Cappe, Robert and Ryden (2003).
 
        \eqn{\mbox{\hspace{6in}}}{\cr}
	If \code{ce}\eqn{=c} is a single value then
	\eqn{c_1=\cdots=c_{K_{max}}=c.}{c[1]=\ldots=c[K[max]]=c.}

        \eqn{\mbox{\hspace{6in}}}{\cr}
	If \code{ce}\eqn{=\boldsymbol{c}}{c} is a vector of length
	\eqn{K_{max}}{K[max]} then the \eqn{j}-th element of \code{ce}
	gives \eqn{c_j}{c[j]} \eqn{(j=1,\dots,K_{max}).}{(j=1,\ldots,K[max]).}      
      }
      \item{D}{a numeric vector or matrix which specifies
        \eqn{\boldsymbol{D}_1, \dots, \boldsymbol{D}_{K_{max}}}{D[1],
	\ldots, D[K[max]]} (prior variances or covariance matrices
        of the mixture means
	\eqn{\boldsymbol{\mu}_1,\dots,\boldsymbol{\mu}_{K_{max}}}{mu[1], \ldots,
          mu[K[max]]} when \code{priormuQ} is
	\dQuote{\code{independentC}}.)
	Its \bold{default} value is a diagonal matrix with squared ranges
	of each column of \code{init$y} on a diagonal.

        \eqn{\mbox{\hspace{6in}}}{\cr}	  
        If \eqn{p=1} and \code{D}\eqn{=d} is a single value
        then \eqn{d_1=\cdots=d_{K_{max}} =
	d.}{d[1]=\ldots=d[K[max]] = d.}
	  
        \eqn{\mbox{\hspace{6in}}}{\cr}
	If \eqn{p=1} and \code{D}\eqn{=\boldsymbol{d}}{=d} is a vector of length
	\eqn{K_{max}}{K[max]} then the \eqn{j}-th element of \code{D}
	gives \eqn{d_j}{d[j]} \eqn{(j=1,\dots,K_{max}).}{(j=1,\ldots,K[max]).}
	  
        \eqn{\mbox{\hspace{6in}}}{\cr}	  
        If \eqn{p>1} and \code{D}\eqn{=\boldsymbol{D}}{=D} is a
	\eqn{p\times p}{p x p} matrix
	then \eqn{\boldsymbol{D}_1=\cdots=\boldsymbol{D}_{K_{max}} =
	  \boldsymbol{D}.}{D[1]=\ldots=D[K[max]] = D.}

        \eqn{\mbox{\hspace{6in}}}{\cr}
	If \eqn{p>1} and \code{D} is a \eqn{(K_{max}\cdot p) \times p}{(K[max]*p)
        x p} matrix then the the first \eqn{p} rows of \code{D}
        give \eqn{\boldsymbol{D}_1}{D[1]}, rows \eqn{p+1,\ldots,2p} of
        \code{D} give \eqn{\boldsymbol{D}_2}{D[2]} etc.
      }
      \item{zeta}{degrees of freedom \eqn{\zeta}{zeta} for the Wishart
	prior on the inverted mixture variances
	\eqn{\boldsymbol{Q}_1,\dots,\boldsymbol{Q}_{K_{max}}.}{Q[1],
	  \ldots,Q[K[max]].}.
	It must be higher then \eqn{p-1}. Its \bold{default} value is
	\eqn{p + 1}.	
      }
      \item{g}{a value or a vector of length \eqn{p} with the shape
	parameters \eqn{g_1,\dots,g_p}{g[1],\ldots,g[p]} for the Gamma hyperpriors on
	\eqn{\gamma_1,\dots,\gamma_p}{gamma[1],\ldots,gamma[p]}.	
	It must be positive. Its \bold{default} value is a vector
	\eqn{(0.2,\dots,0.2)'}{(0.2,\ldots,0.2)'}.
      }
      \item{h}{a value or a vector of length \eqn{p} with the rate
	parameters \eqn{h_1,\dots,h_p}{h[1],\ldots,h[p]}
	for the Gamma hyperpriors on \eqn{\gamma_1,\dots,\gamma_p}{gamma[1],\ldots,gamma[p]}.
	It must be positive. Its \bold{default} value is a vector containing
	\eqn{10/R_l^2}{10/R[l]^2}, where
	\eqn{R_l}{R[l]} is a range of the \eqn{l}-th column of \code{init$y}.
      }	
    }  
  }
  \item{init}{a list with the initial values for the MCMC. All initials
    can be determined by the program if they are not specified. The list
    may have the following components:
    \describe{
      \item{y}{a numeric vector or matrix with the initial values for the latent censored observations.}
      \item{K}{a numeric value with the initial value for the number of mixture components.}
      \item{w}{a numeric vector with the initial values for the mixture
	weights.}
      \item{mu}{a numeric vector or matrix with the initial values for
	the mixture means.}
      \item{Sigma}{a numeric vector or matrix with the initial values
	for the mixture variances.}
      \item{Li}{a numeric vector with the initial values for the Colesky
	decomposition of the mixture inverse variances.}
      \item{gammaInv}{a numeric vector with the initial values for the
	inverted components of the hyperparameter \eqn{\boldsymbol{\gamma}}{gamma}.}
      \item{r}{a numeric vector with the initial values for the mixture
	allocations.}      
    }      
  }
  \item{init2}{a list with the initial values for the second chain
    needed to estimate the penalized expected deviance of Plummer
    (2008). The list \code{init2} has the same structure as the list
    \code{init}. All initials in \code{init2} can be determined by the
    program (differently than the values in \code{init}) if they are not
    specified.

    Ignored when \code{PED} is \code{FALSE}.
  }  
  \item{RJMCMC}{a list with the parameters needed to run reversible jump
    MCMC for mixtures with varying number of components. It does not
    have to be specified if the number of components is fixed. Most of
    the parameters can be determined by the program if they are not
    specified. The list may have the following components:
    \describe{
      \item{Paction}{probabilities (or proportionalit constants) which
	are used to choose an action of the sampler within each
	iteration of MCMC to update the mixture related parameters. Let
	\code{Paction} = \eqn{(p_1,\,p_2,\,p_3)'}{(p[1], p[2],
	  p[3])'}. Then with probability \eqn{p_1}{p[1]} only steps
	assuming fixed \eqn{k} (number of mixture components) are
	performed, with probability \eqn{p_2}{p[2]} split-combine move
	is proposed and with probability \eqn{p_3}{p[3]} birth-death
        move is proposed.

	If not specified (default) then in each iteration of MCMC, all
	sampler actions are performed.
%        Default value \code{Paction} =
%	\eqn{(1/3,\,1/3,\,1/3)'}{(1/3, 1/3, 1/3)'} corresponds to
%	Richardson and Green (1997) where all actions were ued in each iteration.
%        Default value \code{Paction} =
%	\eqn{(0.2,\,0.7,\,0.1)'}{(0.2, 0.7, 0.1)'} follows suggestion of
%	Dellaportas and Papageorgiou (2006).
      }
      \item{Psplit}{a numeric vector of length \code{prior$Kmax} giving
	conditional probabilities of the split move given \eqn{k} as
	opposite to the combine move.

        Default value is
	  \eqn{(1,\,0.5,\ldots,0.5,\,0)'}{(1, 0.5, \ldots, 0.5, 0)'}.
      }
      \item{Pbirth}{a numeric vector of length \code{prior$Kmax} giving
	conditional probabilities of the birth move given \eqn{k} as
	opposite to the death move.
	
        Default value is
	  \eqn{(1,\,0.5,\ldots,0.5,\,0)'}{(1, 0.5, \ldots, 0.5, 0)'}.
      }
      \item{par.u1}{a two component vector with parameters of the beta
	distribution used to generate an auxiliary value
	\eqn{u_1}{u[1]}.

	A default value is \code{par.u1} = \eqn{(2,\,2)'}{(2, 2)'},
	i.e., \eqn{u_1 \sim \mbox{Beta}(2,\,2).}{u[1] ~ Beta(2, 2).}	
      }
      \item{par.u2}{a two component vector (for \eqn{p=1}) or a matrix
	(for \eqn{p > 1})
	with two columns with parameters of the distributions of the
	auxiliary values
	\eqn{u_{2,1},\ldots,u_{2,p}}{u[2,1],\ldots,u[2,p]} in rows.

        A default value leads to
	\eqn{u_{2,d} \sim \mbox{Unif}(-1,\,1)\; (d=1,\ldots,p-1),}{%
	  u[2,d] ~ Unif(-1, 1) (d=1,\ldots,p-1),}
	\eqn{u_{2,p} \sim \mbox{Beta}(1,\,2p).}{u[2,p] ~ Beta(1, 2p).}
	
      }
      \item{par.u3}{a two component vector (for \eqn{p=1}) or a matrix
	(for \eqn{p > 1})
	with two columns with parameters of the distributions of the
	auxiliary values
	\eqn{u_{3,1},\ldots,u_{3,p}}{u[3,1],\ldots,u[3,p]} in rows.

        A default value leads to
	\eqn{u_{3,d} \sim \mbox{Unif}(0,\,1)\; (d=1,\ldots,p-1),}{%
	  u[3,d] ~ Unif(0, 1) (d=1,\ldots,p-1),}
	\eqn{u_{3,p} \sim \mbox{Beta}(1,\,p),}{u[3,p] ~ Beta(1, p).}	
      }      
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
  \item{PED}{a logical value which indicates whether the penalized
    expected deviance (see Plummer, 2008 for more details)
    is to be computed (which requires two parallel
    chains). If not specified, \code{PED} is set to \code{TRUE}
    for models with fixed number of components and is set to
    \code{FALSE} for models with numbers of components estimated using RJ-MCMC.
  }
  \item{keep.chains}{logical. If \code{FALSE}, only summary statistics
    are returned in the resulting object. This might be useful in the
    model searching step to save some memory.}  
  \item{onlyInit}{logical. If \code{TRUE} then the function only
    determines parameters of the prior distribution, initial values,
    values of \code{scale} and
    parameters for the reversible jump MCMC.}
  \item{dens.zero}{a small value used instead of zero when computing
    deviance related quantities.}
  \item{x}{an object of class \code{NMixMCMC} or \code{NMixMCMClist} to
    be printed.}
  \item{dic}{logical which indicates whether DIC should be printed. By
    default, DIC is printed only for models with a fixed number of
    mixture components.
  }
  \item{ped}{logical which indicates whether PED should be printed. By
    default, PED is printed only for models with a fixed number of
    mixture components.
  }
  \item{parallel}{a logical value which indicates whether parallel
    computation (based on a package \code{parallel}) should be used when
    running two chains for the purpose of \code{PED} calculation}
  \item{cltype}{optional argument applicable if \code{parallel} is
    \code{TRUE}. If \code{cltype} is given, it is passed as the
    \code{type} argument into the call to \code{\link[parallel]{makeCluster}}.}  
  \item{\dots}{additional arguments passed to the default \code{print} method.}
}
\value{
  An object of class \code{NMixMCMC} or class \code{NMixMCMClist}.
  Object of class \code{NMixMCMC} is returned if \code{PED} is
  \code{FALSE}. Object of class \code{NMixMCMClist} is returned if
  \code{PED} is \code{TRUE}. 
}
\section{Object of class NMixMCMC}{
  Objects of class \code{NMixMCMC} have the following components:
  \describe{  %%% describe
  \item{iter}{index of the last iteration performed.}
  \item{nMCMC}{used value of the argument \code{nMCMC}.}
  \item{dim}{dimension \eqn{p} of the distribution of data}
  \item{nx_w}{number of levels of a factor covariate on mixture weights
    (equal to 1 if there were no covariates on mixture weights)}
%  \item{x_w}{a factor covariate on mixture weights stored as integers 0,
%	1, \ldots. Only a (scalar) value of 0 if no covariates on
%	mixture weights}
  \item{prior}{a list containing the used value of the argument \code{prior}.}
  \item{init}{a list containing the used initial values for the MCMC
    (the first iteration of the burn-in).
  }
  \item{state.first}{a list having the components labeled
    \code{y}, \code{K}, \code{w}, \code{mu}, \code{Li}, \code{Q}, \code{Sigma},
    \code{gammaInv}, \code{r} containing the values of
    generic parameters at the first stored (after burn-in) iteration of the MCMC.
  }  
  \item{state.last}{a list having the components labeled
    \code{y}, \code{K}, \code{w}, \code{mu}, \code{Li}, \code{Q}, \code{Sigma},
    \code{gammaInv}, \code{r} containing the last sampled values of
    generic parameters.
  }  
  \item{RJMCMC}{a list containing the used value of the argument \code{RJMCMC}.}  
  \item{scale}{a list containing the used value of the argument \code{scale}.}
  \item{freqK}{frequency table of \eqn{K} based on the sampled chain.}
  \item{propK}{posterior distribution of \eqn{K} based on the sampled chain.}
  \item{DIC}{a \code{data.frame} having columns labeled
    \code{DIC}, \code{pD}, \code{D.bar}, \code{D.in.bar} containing
    values used to compute deviance information criterion
    (DIC). Currently only \eqn{DIC_3}{DIC[3]} of Celeux et al. (2006) is
    implemented.
  }
  \item{moves}{a \code{data.frame} which summarizes the acceptance
    probabilities of different move types of the sampler.}
  \item{K}{numeric vector with a chain for \eqn{K} (number of mixture components).}
  \item{w}{numeric vector or matrix with a chain for \eqn{w} (mixture
    weights). It is a matrix with \eqn{K} columns when \eqn{K} is
    fixed. Otherwise, it is a vector with weights put sequentially after
    each other.}
  \item{mu}{numeric vector or matrix with a chain for \eqn{\mu}{mu} (mixture
    means). It is a matrix with \eqn{p\cdot K}{p*K} columns when \eqn{K} is
    fixed. Otherwise, it is a vector with means put sequentially after
    each other.}
  \item{Q}{numeric vector or matrix with a chain for lower triangles of \eqn{\boldsymbol{Q}}{Q} (mixture
    inverse variances). It is a matrix with \eqn{\frac{p(p+1)}{2}\cdot K}{(p*(p+1)2)*K}
    columns when \eqn{K} is fixed. Otherwise, it is a vector with lower
    triangles of \eqn{\boldsymbol{Q}}{Q} matrices put sequentially after each other.}
  \item{Sigma}{numeric vector or matrix with a chain for lower triangles of \eqn{\Sigma}{Sigma} (mixture
    variances). It is a matrix with \eqn{\frac{p(p+1)}{2}\cdot K}{(p*(p+1)2)*K}
    columns when \eqn{K} is fixed. Otherwise, it is a vector with lower
    triangles of \eqn{\Sigma}{Sigma} matrices put sequentially after each other.}
  \item{Li}{numeric vector or matrix with a chain for lower triangles of
    Cholesky decompositions of \eqn{\boldsymbol{Q}}{Q} matrices.
    It is a matrix with \eqn{\frac{p(p+1)}{2}\cdot K}{(p*(p+1)2)*K}
    columns when \eqn{K} is fixed. Otherwise, it is a vector with lower
    triangles put sequentially after each other.}
  \item{gammaInv}{matrix with \eqn{p} columns with a chain for inverses
    of the hyperparameter \eqn{\boldsymbol{\gamma}}{gamma}.
  }
  \item{order}{numeric vector or matrix with order indeces of mixture
    components related to artificial identifiability constraint
    defined by a suitable re-labeling algorithm
    (by default, simple ordering of the first component of the mixture
    means is used).

    It is a matrix with \eqn{K} columns when \eqn{K} is
    fixed. Otherwise it is a vector with orders put sequentially after
    each other.}
  \item{rank}{numeric vector or matrix with rank indeces of mixture
    components. related to artificial identifiability constraint
    defined by a suitable re-labeling algorithm
    (by default, simple ordering of the first component of the mixture
    means is used).

    It is a matrix with \eqn{K} columns when \eqn{K} is
    fixed. Otherwise it is a vector with ranks put sequentially after
    each other.}  
  \item{mixture}{\code{data.frame} with columns labeled
    \code{y.Mean.*}, \code{y.SD.*}, \code{y.Corr.*.*},
    \code{z.Mean.*}, \code{z.SD.*}, \code{z.Corr.*.*} containing the
    chains for the means, standard deviations and correlations of the
    distribution of the original (\code{y}) and scaled (\code{z}) data
    based on a normal mixture at each iteration.
  }
  \item{deviance}{\code{data.frame} with columns labeles
    \code{LogL0}, \code{LogL1}, \code{dev.complete}, \code{dev.observed}
    containing the chains of quantities needed to compute DIC.
  }  
  \item{pm.y}{a \code{data.frame} with \eqn{p} columns with posterior
    means for (latent) values of observed data (useful when there is
    censoring).}
  \item{pm.z}{a \code{data.frame} with \eqn{p} columns with posterior
    means for (latent) values of scaled observed data (useful when there is censoring).}  
  \item{pm.indDev}{a \code{data.frame} with columns labeled
    \code{LogL0}, \code{LogL1}, \code{dev.complete},
    \code{dev.observed}, \code{pred.dens} containing posterior means of
    individual contributions to the deviance.
  }
  \item{pred.dens}{a numeric vector with the predictive density of the
    data based on the MCMC sample evaluated at data points.

    Note that when there is censoring, this is not exactly the
    predictive density as it is computed as the average of densities at
    each iteration evaluated at sampled values of latent observations at
    iterations.
  }
  \item{poster.comp.prob_u}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob_u} is a matrix with \eqn{K} columns and \eqn{n}
    rows with estimated posterior component probabilities
    -- posterior means of the components of the underlying 0/1
    allocation vector.

    \bold{WARNING:} By default,
    the labels of components are based on artificial
    identifiability constraints based on ordering of the mixture means
    in the first margin. Very often, such identifiability constraint is
    not satisfactory!
    
    %These can be used for possible clustering of the subjects.
  }
  \item{poster.comp.prob_b}{a matrix which is present in the output object
    if the number of mixture components in the distribution of random
    effects is fixed and equal to \eqn{K}. In that case,
    \code{poster.comp.prob_b} is a matrix with \eqn{K} columns and \eqn{n}
    rows  with estimated posterior component probabilities
    -- posterior mean over model parameters.

    \bold{WARNING:} By default, the labels of components are based on artificial
    identifiability constraints based on ordering of the mixture means
    in the first margin. Very often, such identifiability constraint is
    not satisfactory!
    
    %These can also be used for possible clustering of the subjects
  }    
  \item{summ.y.Mean}{Posterior summary statistics based on chains stored
    in \code{y.Mean.*} columns of the \code{data.frame} \code{mixture}.}
  \item{summ.y.SDCorr}{Posterior summary statistics based on chains
    stored in \code{y.SD.*} and \code{y.Corr.*.*} columns of the
    \code{data.frame} \code{mixture}.}
  \item{summ.z.Mean}{Posterior summary statistics based on chains stored
    in \code{z.Mean.*} columns of the \code{data.frame} \code{mixture}.}
  \item{summ.z.SDCorr}{Posterior summary statistics based on chains
    stored in \code{z.SD.*} and \code{z.Corr.*.*} columns of the \code{data.frame} \code{mixture}.}  
  \item{poster.mean.w}{a numeric vector with posterior means of mixture
    weights after re-labeling. It is computed only if \eqn{K} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}  
  \item{poster.mean.mu}{a matrix with posterior means of mixture
    means after re-labeling. It is computed only if \eqn{K} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Q}{a list with posterior means of mixture inverse
    variances after re-labeling. It is computed only if \eqn{K} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Sigma}{a list with posterior means of mixture
    variances after re-labeling. It is computed only if \eqn{K} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{poster.mean.Li}{a list with posterior means of Cholesky
    decompositions of mixture inverse variances after re-labeling.
    It is computed only if \eqn{K} is fixed
    and even then I am not convinced that these are useful posterior
    summary statistics (see label switching problem mentioned above).
    In any case, they should be used with care.}
  \item{relabel}{a list which specifies the algorithm used to re-label
    the MCMC output to compute \code{order}, \code{rank}, \code{poster.comp.prob_u},
    \code{poster.comp.prob_b}, \code{poster.mean.w},
    \code{poster.mean.mu}, \code{poster.mean.Q},
    \code{poster.mean.Sigma}, \code{poster.mean.Li}. 
  }
  \item{Cpar}{a list with components useful to call underlying C++ functions (not
    interesting for ordinary users).}
  }
}
\section{Object of class NMixMCMClist}{
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
  \item{popt}{contributions to the unweighted penalty from each observation.}
  \item{wpopt}{contributions to the weighted penalty from each observation.}
  \item{inv.D}{for each observation, number of iterations (in both chains), where the
    deviance was in fact equal to infinity (when the corresponding
    density was lower than \code{dens.zero}) and was not taken into account when
    computing \code{D.expect}.}
  \item{inv.popt}{for each observation, number of iterations, where the
    penalty was in fact equal to infinity and was not taken into account
    when computing \code{popt}.
  }
  \item{inv.wpopt}{for each observation, number of iterations, where the
    importance sampling weight was in fact equal to infinity and was not taken into account
    when computing \code{wpopt}.}
  \item{sumISw}{for each observation, sum of importance sampling
    weights.}
  }
}
\references{
  Celeux, G., Forbes, F., Robert, C. P., and Titterington, D. M. (2006).
  Deviance information criteria for missing data models.
  \emph{Bayesian Analysis}, \bold{1}(4), 651--674.
  
  Cappé, Robert and Rydén (2003).
  Reversible jump, birth-and-death and more general continuous time Markov chain Monte Carlo samplers.
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{65}(3), 679--700.
  
  %Dellaportas, P. and Papageorgiou, I. (2006).
  %Multivariate mixtures of normals with unknown number of components.
  %\emph{Statistics and Computing}, \bold{16}, 57--68.

  Diebolt, J. and Robert, C. P. (1994).
  Estimation of finite mixture distributions through Bayesian sampling.
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{56}(2), 363--375.

  Jasra, A., Holmes, C. C., and Stephens, D. A. (2005).
  Markov chain Monte Carlo methods and the label switching problem in
  Bayesian mixture modelling.
  \emph{Statistical Science}, \bold{20}(1), 50--67.
  
  Komárek, A. (2009).
  A new R package for Bayesian estimation of multivariate normal mixtures allowing for selection 
  of the number of components and interval-censored data.
  \emph{Computational Statistics and Data Analysis}, \bold{53}(12), 3932--3947.

  Plummer, M. (2008).
  Penalized loss functions for Bayesian model comparison.
  \emph{Biostatistics}, \bold{9}(3), 523--539.
  
  Richardson, S. and Green, P. J. (1997).
  On Bayesian analysis of mixtures with unknown number of components
  (with Discussion).
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{59}(4), 731--792.

  Spiegelhalter, D. J.,Best, N. G., Carlin, B. P., and van der Linde, A. (2002).
  Bayesian measures of model complexity and fit (with Discussion).
  \emph{Journal of the Royal Statistical Society, Series B}, \bold{64}(4), 583--639.
}
\seealso{
  \code{\link{NMixPredDensMarg}}, \code{\link{NMixPredDensJoint2}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\examples{
\dontrun{
## See also additional material available in 
## YOUR_R_DIR/library/mixAK/doc/
## or YOUR_R_DIR/site-library/mixAK/doc/
## - files Galaxy.R, Faithful.R, Tandmob.R and
## https://www2.karlin.mff.cuni.cz/~komarek/software/mixAK/Galaxy.pdf
## https://www2.karlin.mff.cuni.cz/~komarek/software/mixAK/Faithful.pdf
## https://www2.karlin.mff.cuni.cz/~komarek/software/mixAK/Tandmob.pdf
##

## ==============================================

## Simple analysis of Anderson's iris data
## ==============================================
library("colorspace")

data(iris, package="datasets")
summary(iris)
VARS <- names(iris)[1:4]
#COLS <- rainbow_hcl(3, start = 60, end = 240)
COLS <- c("red", "darkblue", "darkgreen")
names(COLS) <- levels(iris[, "Species"])

### Prior distribution and the length of MCMC
Prior <- list(priorK = "fixed", Kmax = 3)
nMCMC <- c(burn=5000, keep=10000, thin=5, info=1000)

### Run MCMC
set.seed(20091230)
fit <- NMixMCMC(y0 = iris[, VARS], prior = Prior, nMCMC = nMCMC)

### Basic posterior summary
print(fit)

### Univariate marginal posterior predictive densities
### based on chain #1
pdens1 <- NMixPredDensMarg(fit[[1]], lgrid=150)
plot(pdens1)
plot(pdens1, main=VARS, xlab=VARS)

### Bivariate (for each pair of margins) predictive densities
### based on chain #1
pdens2a <- NMixPredDensJoint2(fit[[1]])
plot(pdens2a)

plot(pdens2a, xylab=VARS)
plot(pdens2a, xylab=VARS, contour=TRUE)

### Determine the grid to compute bivariate densities
grid <- list(Sepal.Length=seq(3.5, 8.5, length=75),
             Sepal.Width=seq(1.8, 4.5, length=75),
             Petal.Length=seq(0, 7, length=75),
             Petal.Width=seq(-0.2, 3, length=75))
pdens2b <- NMixPredDensJoint2(fit[[1]], grid=grid)
plot(pdens2b, xylab=VARS)

### Plot with contours
ICOL <- rev(heat_hcl(20, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)))
oldPar <- par(mfrow=c(2, 3), bty="n")
for (i in 1:3){
  for (j in (i+1):4){
    NAME <- paste(i, "-", j, sep="")
    MAIN <- paste(VARS[i], "x", VARS[j])
    image(pdens2b$x[[i]], pdens2b$x[[j]], pdens2b$dens[[NAME]], col=ICOL,
          xlab=VARS[i], ylab=VARS[j], main=MAIN)
    contour(pdens2b$x[[i]], pdens2b$x[[j]], pdens2b$dens[[NAME]], add=TRUE, col="brown4")
  }  
}  

### Plot with data
for (i in 1:3){
  for (j in (i+1):4){
    NAME <- paste(i, "-", j, sep="")
    MAIN <- paste(VARS[i], "x", VARS[j])
    image(pdens2b$x[[i]], pdens2b$x[[j]], pdens2b$dens[[NAME]], col=ICOL,
          xlab=VARS[i], ylab=VARS[j], main=MAIN)
    for (spec in levels(iris[, "Species"])){
      Data <- subset(iris, Species==spec)
      points(Data[,i], Data[,j], pch=16, col=COLS[spec])
    }  
  }  
}  

### Set the graphical parameters back to their original values
par(oldPar)

### Clustering based on posterior summary statistics of component allocations
### or on the posterior distribution of component allocations
### (these are two equivalent estimators of probabilities of belonging
###  to each mixture components for each observation)
p1 <- fit[[1]]$poster.comp.prob_u
p2 <- fit[[1]]$poster.comp.prob_b

### Clustering based on posterior summary statistics of mixture weight, means, variances
p3 <- NMixPlugDA(fit[[1]], iris[, VARS])
p3 <- p3[, paste("prob", 1:3, sep="")]

  ### Observations from "setosa" species (all would be allocated in component 1)
apply(p1[1:50,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p2[1:50,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p3[1:50,], 2, quantile, prob=seq(0, 1, by=0.1))

  ### Observations from "versicolor" species (almost all would be allocated in component 2)
apply(p1[51:100,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p2[51:100,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p3[51:100,], 2, quantile, prob=seq(0, 1, by=0.1))

  ### Observations from "virginica" species (all would be allocated in component 3)
apply(p1[101:150,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p2[101:150,], 2, quantile, prob=seq(0, 1, by=0.1))
apply(p3[101:150,], 2, quantile, prob=seq(0, 1, by=0.1))
}
}
\keyword{multivariate}
\keyword{survival}
