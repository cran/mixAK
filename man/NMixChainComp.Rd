\name{NMixChainComp}
\alias{NMixChainComp}
\alias{NMixChainComp.default}
\alias{NMixChainComp.NMixMCMC}
\alias{NMixChainComp.GLMM_MCMC}
\title{
  Chains for mixture parameters
}
\description{
  This function returns chains for parameters derived from the
  (re-labeled) mixture weights, means, covariance matrices.

  First, mixture means and shifted-scaled to the original (data) scale, mixture
  covariance matrices are scaled to the original (data) scale
  (see argument \code{scale} in \code{\link{NMixMCMC}} function
  or argument \code{scale.b} in \code{\link{GLMM_MCMC}}).
  Possible derived parameters are standard deviations and correlation coefficients.
}
\usage{
NMixChainComp(x, relabel = TRUE, param)

\method{NMixChainComp}{default}(x, relabel = TRUE, param)

\method{NMixChainComp}{NMixMCMC}(x, relabel = TRUE,
   param = c("w", "mu", "var", "sd", "cor", "Sigma", "Q", "Li"))

\method{NMixChainComp}{GLMM_MCMC}(x, relabel = TRUE,
   param = c("w_b", "mu_b", "var_b", "sd_b", "cor_b", "Sigma_b", "Q_b", "Li_b"))
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} or \code{GLMM_MCMC}.}
  \item{relabel}{a logical argument indicating whether the chains are to
    be returned with components being re-labeled (see
    \code{\link{NMixRelabel}}) or whether the chains are to be returned as
    originally sampled.}
  \item{param}{a character string indicating which sample is to be
    returned:
    \describe{
      \item{w, w_b}{mixture weights;}
      \item{mu, mu_b}{mixture means;}
      \item{var, var_b}{mixture variances;}
      \item{sd, sd_b}{mixture standard deviations;}
      \item{cor, cor_b}{correlations derived from the mixture covariance
        matrices;}
      \item{Sigma, Sigma_b}{mixture covariance matrices (their lower triangles);}
      \item{Q, Q_b}{mixture inverted covariance matrices (their lower triangles);}
      \item{Li, Li_b}{Cholesky factors (their lower triangles) of the
	mixture inverted covariance matrices.}      
    }  
  }
}
\value{
  A matrix with sampled values in rows, parameters in columns.
}
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}  
\keyword{multivariate}
\keyword{smooth}
