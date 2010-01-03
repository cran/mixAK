\name{GLMM_MCMCifit}
\alias{GLMM_MCMCifit}
\title{
  Initial (RE)ML fits for the GLMM_MCMC function
}
\description{
  This is a help function for \code{\link{GLMM_MCMC}} function.
  Besides initial (RE)ML fits, the function created variables
  derived from the design matrices.

  THIS FUNCTION IS NOT TO BE CALLED BY ORDINARY USERS.
}  
\usage{
GLMM_MCMCifit(do.init, na.complete,
    y, dist, id, time, x, z, random.intercept,
    xempty, zempty, Rc, Rd, p, p_fi, q, q_ri, lbeta, dimb)
}
\arguments{
  \item{do.init}{logical value indicating whether initial (RE)ML fits
    should be done}
  \item{na.complete}{logical value. If \code{TRUE} then
    the function removes rows containing NA's from \code{y}, \code{id}, \code{x},
    \code{z} whenever there is at least one missing value for arbitrary
    response. If \code{FALSE} then the missing values are removed
    response by response, i.e., different response variables may have
    different numbers of observations.}
  \item{y}{see output element \code{y} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{dist}{see argument\code{dist} of \code{\link{GLMM_MCMC}}
    function}
  \item{id}{see output element \code{id} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{time}{see argument \code{time} of \code{GLMM_longitClust}}
  \item{x}{see output element \code{x} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{z}{see output element \code{z} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{random.intercept}{see output element \code{random.intercept} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{xempty}{see output element \code{xempty} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{zempty}{see output element \code{zempty} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{Rc}{see output element \code{Rc} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{Rd}{see output element \code{Rd} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{p}{see output element \code{p} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{p_fi}{see output element \code{p_fi} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{q}{see output element \code{q} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{q_ri}{see output element \code{q_ri} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{lbeta}{see output element \code{lbeta} of \code{\link{GLMM_MCMCdata}}
    function}
  \item{dimb}{see output element \code{dimb} of \code{\link{GLMM_MCMCdata}}
    function}  
}
\value{
  A list with the following components (some of them not included if
  \code{do.init} is \code{FALSE}):
  \item{Y}{a list of length R with observations really used in fitting process
     (after removal of missing values)}
  \item{ID}{a list of length R with id's corresponding to Y}
  \item{time}{a vector \code{time} upon removal of missing values} 
  \item{x}{a list resulting from the original argument x after removal of
     observations with some missing information
     additionaly, intercept column is added if fixed intercept included in the model}
  \item{z}{a list resulting from the original argument z after removal of
     observations with some missing information
     additionaly, intercept column is added if random intercept included in the model}
  \item{I}{number of clusters in the original data (before removing NA's)}
  \item{n}{a list of length R, each component is a vector or length I (may contain zeros
     if some cluster disappears for particular response due to NA's)}
  \item{Cn}{vectorized n}
  \item{sumCn}{sum(Cn) = total number of observations}
  \item{Cy_c}{vector with continuous response to be passed to C++, equal to 0 if there is no continuous response}
  \item{Cy_d}{vector with discrete response to be passed to C++, equal to 0 if there is no discrete response}
  \item{CX}{vector containing X matrices (without ones for possible intercept)
     to be passed to C++, equal to 0 if there are no X matrices}
  \item{CZ}{vector containing Z matrices (without ones for possible intercept)
     to be passed to C++, equal to 0 if there are no Z matrices}
%  \item{CXtX}{vector containing lower triangles of X[s]'X[s] for each response to be passed to C++,
%     where X[s] contains also column of ones if there is a fixed intercept in the model,
%     equal to 0 if there are no beta coefficients in the model}
%  \item{CZitZi}{vector containig lower triangles of matrices Z[i,s]'Z[i,s] for each cluster and each response
%     to be passed to C++, where Z[i,s] contains also column of ones if there is a random intercept in the model,
%     Remark:   for cluster without any observations of response s, Z[i,s] is considered
%               to be a matrix of zeros
%               => for CZitZi, matrices Z[i,s]'Z[i,s] are for all i of the same dimension 
%     ordering: cluster 1: ZitZi matrices for s=1,...,R, cluster 2: ZitZi matrices for s=1,...,R, etc.,
%     equal to 0 if there are no random effects in the model}
  \item{iintcpt}{data.frame(Est, SE) with estimated intercepts and their SE, R rows,
      row equal to (0, 0) if there is no fixed intercept for particular response}
  \item{ifixef}{a list of length R, each component is equal to 0 if there are no fixed effects for particular response,
      and is equal to data.frame(Est, SE) if there are fixed effects}
  \item{isigma}{vector of length R, equal to 0 for discrete response, equal to estimated residual
      standard deviation for continuous response}
  \item{iEranef}{a list of length R, each component is equal to 0 if there are no random effects for particular response,
      and is equal to data.frame(Est, SE) with estimated means of the random effects
      and their std. errors if there are random effects}
  \item{iSDranef}{a list of length R, each component is equal to 0 if there are no random effects for particular response,
      and is equal to a vector with estimated standard deviations of the random effects if there are random effects}
  \item{ib}{a list of length R, each component is equal to 0 if there are no random effects for particular response,
      and a matrix with EB estimates of random effects shifted by their estimated mean if there are random effects}
  \item{is.intcpt}{logical vector of length R}
  \item{is.fixef}{logical vector of length R}
  \item{is.ranef}{logical vector of length R}  
  \item{is.sigma}{logical vector of length R}
  \item{ibMat}{matrix with initial values of random effects (EB estimates from (RE)ML fits)}
  \item{iEranefVec}{vector with estimated means of random effects}
  \item{iSEranefVec}{vector with standard errors of estimated means of random effects}
  \item{iSDranefVec}{vector with estimated standard deviations of random effects}
  \item{ibeta}{vector with initial values of beta's (including fixed intercepts)}
}
\seealso{
  \code{\link{GLMM_MCMC}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{internal}
