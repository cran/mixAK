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

  \method{tracePlots}{NMixMCMC}(x, param=c("Emix", "SDmix", "Cormix", "w", "mu", "gammaInv"),
    auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", \dots)

  \method{tracePlots}{GLMM_MCMC}(x, param=c("beta", "Eb", "SDb", "Corb", "sigma_eps",
    "w_b", "mu_b", "gammaInv_b", "gammaInv_eps"),
    auto.layout=TRUE, xlab="Iteration", ylab, col="slateblue", main="", \dots)
}
\arguments{
  \item{x}{an object of class \code{NMixMCMC} or \code{GLMM_MCMC}.}
  \item{param}{a character string which specifies which sort of
    parameters is to be plotted.}
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
  \code{\link{NMixMCMC}}, \code{\link{GLMM_MCMC}}, \code{\link[coda]{traceplot}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{dplot}
