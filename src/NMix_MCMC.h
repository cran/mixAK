//
//  PURPOSE:   Main functions to run RJ-MCMC to estimate
//             a multivariate density
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/10/2007
//
//  FUNCTIONS:  
//     * NMix_MCMC  29/10/2007:  Start working on it
//                  21/12/2007:  Working version for fixed K, seems to be bugs free
//                  29/01/2008:  Version for varying K with p = 1 seems to work
//                               Version for varying K with p > 1 finished but does not seem to work ;-(
//                  30/12/2009:  Computation of quantities needed for clustering added
//                  27/03/2015:  Mild modification to allow for covariates on mixture weights
//                               (if number of components is fixed)
//                               - argument nxw_xw added
//
// ======================================================================
//
#ifndef _NORMAL_MIXTURE_MCMC_H_
#define _NORMAL_MIXTURE_MCMC_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
#include <R_ext/RS.h>

#include "AK_Basic.h"
#include "AK_Utils.h"
#include "AK_BLAS.h"

#include "Dist_MVN.h"
//#include "Dist_Wishart.h"

#include "NMix.h"
#include "NMix_Utils.h"
#include "NMix_orderComp.h"
#include "NMix_Deviance.h"
#include "NMix_update_sum_Ir_and_sum_Pr_y.h"

#include "NMix_updateAlloc.h"
#include "NMix_updateWeights.h"
#include "NMix_updateMeansVars.h"
#include "NMix_updateHyperVars.h"
#include "NMix_updateCensObs.h"

#include "NMix_RJMCMC_aux_vector_u.h"
#include "NMix_RJMCMCsplit.h"
#include "NMix_RJMCMCcombine.h"
#include "NMix_RJMCMCbirth.h"
#include "NMix_RJMCMCdeath.h"

#include "NMix_PosterMeanMixParam.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_MCMC                                                                                 *****/
/***** ***************************************************************************************** *****/
//
// y0[p, n]        data (matrix stored in column major order)
//                 * exactly observed data
//                 * right-censored observations                              
//                 * left-censored observations
//                 * lower limits of interval-censored observations
//                  
// y1[p, n]        data (matrix stored in column major order)
//                 * upper limits of interval-censored observations
//
// censor[p, n]    censoring indicator for each observation (matrix stored in column major order)
//                 * 0 = right-censored
//                 * 1 = observed
//                 * 2 = left-censored
//                 * 3 = interval-censored
//
// nxw_xw[1 + n]   information on a factor covariate on mixture weights (ADDED ON 20150327)
//                 nxw_xw[0]   = nxw = number of covariate levels
//                 nxw_xw[1:n] = covariate values (numbered 0, 1, ..., nxw - 1)
//
// dimy[2]         dimension of the data
//                 * dimy[0]  = p = dimension of the response
//                 * dimy[+1] = n = number of observations
//
// shiftScale[p+p] shift[p] and scale[p] which was used to get y0 and y1 from the original data
//                 * it is used to calculate chCorr
//
// nMCMC[4]        parameters of the MCMC
//                 * nMCMC[0]  = Mburn = number of burn-in scans (after thinning has been applied)
//                               * can be equal to zero (no burn-in)
//                 * nMCMC[+1] = Mkeep = number of additional scans (after thinning has been applied)
//                                      (scans which are stored for the inference)
//                 * nMCMC[+1] = thinning interval
//                 * nMCMC[+1] = info interval (how often we write an information concerning the MCMC on the screen)
//
// priorInt[3]     integer prior hyperparameters
//                 * priorInt[0]  = type of prior for K (see enum _NMix_type_priorK)
//                                  * 0 -> K is fixed
//                                  * 1 -> uniform on 1, ..., Kmax
//                                  * 2 -> Poisson(lambda) truncated on Kmax
//                 * priorInt[+1] = type of prior for (mu[j], Sigma[j]^{-1}) (see enum _NMix_type_priormuSigma)
//                                  * 0 -> natural conjugate
//                                  * 1 -> independent conjugate
//                 * priorInt[+1] = Kmax
//
// priorDouble[?]  double prior hyperparameters
//                 * priorDouble[0]  = lambda 
//                                     * expectation of the Poisson distribution if the truncated Poisson is used as a prior for K
//                 * priorDouble[+1] = delta 
//                                     * prior "sample size" for the Dirichlet prior on weights
//                 * priorDouble[+p*Kmax] = xi[j] (j=0, Kmax-1) 
//                                          * prior means for the mixture means
//                 * priorDouble[+Kmax]   = c[j] (j=0, Kmax-1)
//                                          * prior scale parameters for the mixture means if natural conjugate prior is used
//                 * priorDouble[+LT(p)*Kmax]  = D[j]^{-1} (j=0, ..., Kmax-1)                 
//                                          * prior inverse variances for the mixture means if independent conjugate prior is used
//                 * priorDouble[+1]      = zeta
//                                          * prior degrees of freedom of the Wishart prior on Sigma[j]^{-1}
//                 * priorDouble[+p]      = g
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1}
//                 * priorDouble[+p]      = h
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1}
//
// priorRJMCMC[?]  parameters related to the RJ-MCMC
//                 * priorRJMCMC[3]    = Paction
//                                       * probabilities to choose an action of the sampler when updating mixture parameters
//                                         with varying K
//                                       * ignored if priorRJMCMCint[0] = 1
//                 * priorRJMCMC[+Kmax] = Psplit
//                                        * probabilities of the split move (compared to the combine move) given k
//                 * priorRJMCMC[+Kmax] = Pbirth
//                                        * probabilities of the birth move (compared to the death move) given k
//                 * priorRJMCMC[+2] = par.u1
//                                     * parameters of the distribution of the auxiliary number u1
//                 * priorRJMCMC[+2*p] = par.u2
//                                       * parameters of the distribution of the auxiliary vector u2
//                 * priorRJMCMC[+2*p] = par.u3
//                                       * parameters of the distribution of the auxiliary vector u3
//
// priorRJMCMCint[1]  integer parameters related to RJ-MCMC
//
//                 * priorRJMCMCint[0]   = actionDeterministic
//                                         if <> 0 then within RJ-MCMC, all sampler actions are used in each iteration
//
// y[p, n]             INPUT: initial values for (censored) observations
//                            * must be equal to the observed values if there is no censoring
//                     OUTPUT:  last sampled value
//
// y_first[p, n]       INPUT:   whatsever
//                     OUTPUT:  value of y corresponding to the first kept (not burn-in, after thinning loop) MCMC iteration
//
// K[1]                INPUT: initial number of mixture components
//                     OUTPUT:  last sampled value
//
// w[Kmax]             INPUT: initial values of the weights
//                     OUTPUT:  last sampled value
//
// mu[p, Kmax]         INPUT: initial values of the mixture means
//                     OUTPUT:  last sampled value
//
// Q[LT(p), Kmax]      INPUT:  whatsever (initial Q is computed from initial Li)
//                     OUTPUT: last sampled value of mixture inverse variances
//
// Sigma[LT(p), Kmax]   INPUT:  whatsever (initial Sigma is computed from initial Li)
//                     OUTPUT: last sampled value of mixture variances
//
// Li[LT(p), Kmax]     INPUT: initial values of the Cholesky decompositions of mixture inverse variances, i.e., Q[j] = Sigma[j]^{-1} = Li[j] %*% t(Li[j])
//                     OUTPUT:  last sampled value
//
// gammaInv[p]         INPUT:  initial values for the inverted diagonal of the matrix Xi (hyperprior above Sigma[j]^{-1})
//                     OUTPUT:  last sampled value
//
// r[n]                INPUT: initial values for allocations (0,...,K-1)
//                     OUTPUT:  last sampled value
//
// r_first[n]          INPUT:   whatsever
//                     OUTPUT:  value of r corresponding to the first kept (not burn-in, after thinning loop) MCMC iteration
//
// chK[M]              INPUT:  whatsever
//                     OUTPUT: sampled values of K
//
// chw[Kmax*M]         INPUT:  whatsever
//                     OUTPUT: sampled values of w (on the first sum(K) places)
//
// chmu[p*Kmax*M]      INPUT:  whatsever
//                     OUTPUT: sampled values of mu's (on the first p*sum(K) places)
//
// chQ[LT(p)*Kmax*M]   INPUT:  whatsever
//                     OUTPUT: sampled values of Q's, where Q[j] = Sigma[j]^{-1}
//                             (on the first LT(p)*sum(K) places)
//
// chSigma[LT(p)*Kmax*M]   INPUT:  whatsever
//                        OUTPUT: sampled values of Sigma's (mixture variances)
//                                (on the first LT(p)*sum(K) places)
//
// chLi[LT(p)*Kmax*M]  INPUT:  whatsever
//                     OUTPUT: sampled values of Li's, where Li[j] is the Cholesky decomposition of Q[j] = Sigma[j]^{-1}
//                             (on the first LT(p)*sum(K) places)
//
// chgammaInv[p, M]    INPUT:  whatsever
//                     OUTPUT: sampled values of inverted gamma's
//
// chorder[Kmax*M]     INPUT:  whatsever
//                     OUTPUT: ordering of the mixture components for each sampled mixture
//                             (on the first sum(K) places)
//                             - C indeces from 0, ..., K-1
//
// chrank[Kmax*M]      INPUT:  whatsever
//                     OUTPUT: ranks of the mixture components for each sampled mixture
//                             (on the first sum(K) places)
//                             - C indeces from 0, ..., K-1
//
// chMean[p, M]        INPUT:  whatsever
//                     OUTPUT: values of the mixture (overall) mean at each iteration
//
// chCorr[LT(p), M]    INPUT:  whatsever
//                     OUTPUT: values of the mixture (overall) standard deviations (diagonal) and correlations (off-diagonal) at each iteration 
//
// chMeanData[p, M]    INPUT:  whatsever
//                     OUTPUT:  values of the mixture (overall) mean at each iteration after taking into account shift and scale,
//                              i.e., this applies to the original data
//
// chCorrData[LT(p), M]    INPUT:  whatsever
//                        OUTPUT:  values of the mixture (overall) standard deviations (diagonal) and correlations (off-diagonal) at each iteration 
//                                 after taking into account shift and scale, i.e., this applies to the original data
//
// chLogL0[M]           INPUT:   whatsever
//                     OUTPUT:   values of the complete data log-likelihood, part 1, at each iteration, i.e., values of
//                               sum_{i=1}^n log(phi(y_i | mu_{r_i}, Sigma_{r_i})))
//
// chLogL1[M]           INPUT:   whatsever
//                     OUTPUT:   values of the complete data log-likelihood, part 2, at each iteration, i.e., values of
//                               sum_{i=1}^n log(w_{r_i})
//
// chDevCompl[M]        INPUT:   whatsever
//                     OUTPUT:   values of the complete data deviance at each iteration 
//                               (eq. (7) in Celeux, Forbes, Robert, Titterington, 2006)
//
// chDevObs[M]          INPUT:   whatsever
//                     OUTPUT:   values of the observed data deviance at each iteraion, i.e., values of
//                               -2 * sum_{i=1}^n log(sum_{j=1}^K w_j * phi(y_i | mu_j, Sigma_j))
//
// chDevCompl_inHat[M]  INPUT:   whatsever
//                     OUTPUT:   values of the quantities needed to compute the second part of DIC_4 (see Celeux, Forbes, Robert, Titterington, 2006)
//                               = -2 * sum_{i=1}^n (log(E[w_{r_i}|...]) + log(y_i | E[mu_{r_i}|...], (E[Q_{r_i} | ...])^{-1})) 
//
// pm_y[p, n]           INPUT:  whatsever
//                     OUTPUT:  posterior mean (based on MCMC) of y
//
// pm_indLogL0[n]    INPUT:  whatsever
//                  OUTPUT:  posterior mean of individual components of the complete data log-likelihood, part 1
//
// pm_indLogL1[n]        INPUT:  whatsever
//                      OUTPUT:  posterior mean of individual components of the complete data log-likelihood, part 2
//
// pm_indDevCompl[n]    INPUT:  whatsever
//                     OUTPUT:  posterior mean of individual components of the complete data deviance
//
// pm_indDevObs[n]      INPUT:  whatsever
//                     OUTPUT:  posterior mean of individual components of the observed data deviance
//
// pm_indDevCompl_inHat[n]      INPUT:  whatsever
//                             OUTPUT:  posterior mean of individual components of the complete data deviance evaluated in "theta{hat}"
//
// pm_pred_dens[n]      INPUT:  whatsever
//                     OUTPUT:  if there is no censoring:  this will be a value of the predictive density evaluated at observed points
//                              if censoring: pred_dens[i] will be a value of (1/M) sum_{l=1}^M f^{(l)}(y_i^{(l)}), where
//                                            y_i^{(l)} is a sampled value of "true" y_i at l-th iteration
//                                            and f^{(l)} is the mixture at the l-th iteration 
//
// pm_w[Kmax]           INPUT:  whatsever
//                     OUTPUT:  posterior mean of mixture weights, computed only if K is FIXED
//                              * before the posterior mean is computed, components are re-ordered to satisfy some constraint
//
// pm_mu[p, Kmax]       INPUT:  whatsever
//                     OUTPUT:  posterior mean of mixture means, computed only if K is FIXED
//                              * before the posterior mean is computed, components are re-ordered to satisfy some constraint
//
// pm_Q[LT(p), Kmax]    INPUT:   whatsever
//                     OUTPUT:  posterior mean of mixture inverse variances, computed only if K is FIXED
//                              * before the posterior mean is computed, components are re-ordered to satisfy some constraint
//
// pm_Sigma[LT(p), Kmax]   INPUT:   whatsever
//                        OUTPUT:  posterior mean of mixture variances, computed only if K is FIXED
//                                 * before the posterior mean is computed, components are re-ordered to satisfy some constraint
//
// pm_Li[LT(p), Kmax]       INPUT:  whatsever
//                         OUTPUT:  posterior mean of Cholesky decompositions of mixture inverse variances, computed only if K is FIXED
//                                  * before the posterior mean is computed, components are re-ordered to satisfy some constraint
//
// sum_Ir[n, Kmax]:    for each observation and each mixture component: sum(r[i] = k),
//                     i = 0, ..., n-1, j = 0, ..., K - 1
//                     from the main part of MCMC (burn-in not included) 
//                     * COMPUTED ONLY WHEN K is FIXED
//                     * COMPONENTS ARE (INTERNALLY) RE-LABELED BEFORE sum(r[i] = k) is computed
//
// sum_Pr_y[n, Kmax]:    for each observation and each mixture component: sum(P(r[i] = k | theta, y))
//                       i = 0, ..., n-1, j = 0, ..., K - 1
//                       from the main part of MCMC (burn-in not included) 
//                       * COMPUTED ONLY WHEN K is FIXED
//                       * COMPONENTS ARE (INTERNALLY) RE-LABELED BEFORE sum(P(r[i] = k | theta, y)) is computed
//
// iter[1]              INPUT:  index of the iteration corresponding to the initial values (usually zero)
//                     OUTPUT:  index of the last performed iteration
//
// nMoveAccept[9]       INPUT:  whatsever
//                     OUTPUT:  
//                         nMoveAccept[0] = nGibbs_K
//                         nMoveAccept[1] = nSplit
//                         nMoveAccept[2] = nCombine
//                         nMoveAccept[3] = nBirth
//                         nMoveAccept[4] = nDeath
//                         nMoveAccept[5] = nAcceptSplit
//                         nMoveAccept[6] = nAcceptCombine
//                         nMoveAccept[7] = nAcceptBirth
//                         nMoveAccept[8] = nAcceptDeath
//                         --> related to the main MCMC only (it does not count the burn-in)
//
// err[1]
//
void
NMix_MCMC(const double* y0,  
          const double* y1,     
          const int* censor,    
	  const int* nxw_xw,      
          const int* dimy,            
          const double* shiftScale,
          const int* nMCMC,  
          const int* priorInt,  
          const double* priorDouble,  
          const double* priorRJMCMC,  
          const int* priorRJMCMCint,
          double* y,  
          double* y_first,          
          int* K,               
          double* w,                 
          double* mu,                
          double* Q,            
          double* Sigma,        
          double* Li,                
          double* gammaInv,      
          int* r,
          int* r_first,
          int* chK,             
          double* chw,          
          double* chmu,           
          double* chQ,          
          double* chSigma,      
          double* chLi,              
          double* chgammaInv,    
          int* chorder,                  
          int* chrank,
          double* chMean,       
          double* chCorr,       
          double* chMeanData,        
          double* chCorrData,
          double* chLogL0,      
          double* chLogL1,      
          double* chDevCompl,        
          double* chDevObs,        
          double* chDevCompl_inHat,  
          double* pm_y,         
          double* pm_indLogL0,  
          double* pm_indLogL1,  
          double* pm_indDevCompl,    
          double* pm_indDevObs,  
          double* pm_indDevCompl_inHat,  
          double* pm_pred_dens,
          double* pm_w,         
          double* pm_mu,          
          double* pm_Q,              
          double* pm_Sigma,      
          double* pm_Li,
          int*    sum_Ir,
          double* sum_Pr_y,
          int* iter,            
          int* nMoveAccept,     
          int* err);

#ifdef __cplusplus
}
#endif

#endif
