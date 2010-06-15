\name{tracePlots}
\alias{tracePlots}
\alias{tracePlots.default}
\alias{tracePlots.NMixMCMC}
\alias{tracePlots.GLMM_MCMC}
\title{
  Traceplots for selected parameters
}
\description{
  This function draws traceplots of selected parameters from the MCMC
  simulations ran using \code{\link{NMixMCMC}} or
  \code{\link{GLMM_MCMC}} functions.
}
\usage{
  tracePlots(x, \dots)

  \method{tracePlots}{default}(x, \dots)

  \method{tracePlots}{NMixMCMC}(x, param=c("Emix", "SDmix", "Cormix", "K", "w", "mu", "sd", "gammaInv"),
    relabel=FALSE, order,
    auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", \dots)

  \method{tracePlots}{GLMM_MCMC}(x, param=c("Deviance", "Cond.Deviance",
    "beta", "Eb", "SDb", "Corb", "sigma_eps",
    "w_b", "mu_b", "sd_b", "gammaInv_b", "gammaInv_eps"),
    relabel=FALSE, order,
    auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} or \code{GLMM_MCMC}.}
  \item{param}{a character string which specifies which sort of
    parameters is to be plotted.
   
    \describe{
      \item{Emix}{overall means (for each margin) of the normal mixture;}
      \item{SDmix}{overall standard deviations (for each margin) of the
	normal mixture;}
      \item{Cormix}{overall correlations (each pair) of the normal
	mixture;}
      \item{K}{number of mixture components;}
      \item{w, w_b}{weights of each of mixture components.
        If \code{relabel} is \code{FALSE}, weights are not re-labeled before plotting;}
      \item{mu, mu_b}{component means (each margin, each mixture component) of the
	normal mixture. The mixture means are shifted and scaled using
	\code{x$scale$shift} and \code{x$scale$scale} before plotting.
	If \code{relabel} is \code{FALSE}, means are not re-labeled before plotting;}
      \item{sd, sd_b}{component standard deviations (each margin, each mixture
	component) of the normal mixture. The mixture standard deviations are scaled using
	\code{x$scale$scale} before plotting.
	If \code{relabel} is \code{FALSE}, standard deviations are
	not re-labeled before plotting;}
      \item{gammaInv, gammaInv_b, gammaInv_eps}{variance
	hyperparameters;}
      \item{Deviance}{deviance (marginal, with random effects integrated
        out) of the GLMM;}
      \item{Cond.Deviance}{conditional deviance (given random effects)
	of the GLMM;}
      \item{beta}{fixed effects of the fitted GLMM;}
      \item{Eb}{overall means (for each margin) of the random effects of
        the fitted GLMM;}
      \item{SDb}{overall standard deviations (for each margin) of the random effects of
        the fitted GLMM;}
      \item{Corb}{overall correlations (each pair) of the distribution
	of the random effects of the fitted GLMM.}
      \item{sigma_eps}{standard deviations of the error terms in the
	(mixed) models for continuous responses.}
    }  
  }
  \item{relabel}{logical value. It indicates whether the chains with
    \code{param} being \code{w}, \code{mu}, \code{sd}, \code{w_b},
    \code{mu_b}, \code{sd_b} should be re-labeled before
    plotting. Re-labelling is given by argument \code{order}. If
    \code{order} is missing then \code{x$order} or \code{x$order_b}
    determines re-labelling.
  }
  \item{order}{a matrix with \eqn{K} columns and \eqn{M} rows where
    \eqn{M} is the length of MCMC. Each row of \code{order} must be a
    permutation of \eqn{(1, \ldots, K)}{(1, ..., K)}.
  }
  \item{auto.layout}{logical value. If \code{TRUE}, the plotting region
    is automatically divided to produce traceplots of all
    parameters. Note that layout must be set up automatically if there
    are more than 28 parameters to be plotted (often the case for
    correlations with \code{param} being \code{Corb} or for mixture
    means with \code{param} being \code{mu_b}).}
  \item{xlab, ylab, col, main}{arguments passed to \code{plot}
    function. They all can be of length one (the value is used on all
    plots) or of length equal to the number of parameters to be
    plotted.}
  \item{\dots}{other arguments passed to \code{plot} function.}
}
\value{
  \code{invisible(x)}
}
\seealso{
  \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}, \code{\link{NMixRelabel}}, \code{\link[coda]{traceplot}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{dplot}
