//
//  PURPOSE:   Main functions to run MCMC to estimate a (multivariate) GLMM
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/07/2009
//
//  FUNCTIONS:  
//     * GLMM_MCMC  03/07/2009:  Start working on it
//                  03/08/2009:  Version allowing for gaussian responses
//                  10/11/2009:  Version allowing for mixed gaussian and discrete response
//                           
// ============================================================================================
//
#ifndef _GLMM_MCMC_H_
#define _GLMM_MCMC_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"
#include "AK_Utils.h"
#include "AK_BSTAT.h"

#include "NMix.h"
#include "NMix_Utils.h"
#include "NMix_orderComp.h"
#include "NMix_Deviance.h"

#include "NMix_updateAlloc.h"
#include "NMix_updateWeights.h"
#include "NMix_updateMeansVars.h"
#include "NMix_updateHyperVars.h"

#include "NMix_PosterMeanMixParam.h"

#include "GLMM.h"
#include "GLMM_Deviance.h"

#include "GLMM_linear_predictors.h"
//#include "GLMM_scale_ZitZi.h"
#include "GLMM_create_SZitZiS.h"
#include "GLMM_create_XtX.h"
#include "GLMM_create_ZS.h"
#include "GLMM_dY_meanY.h"

#include "GLMM_updateVars_eps.h"
#include "GLMM_updateHyperVars_eps.h"
#include "GLMM_updateFixEf.h"
#include "GLMM_updateRanEf.h"
#include "GLMM_updateRanEf_QR.h"


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_MCMC                                                                                 *****/
/***** ***************************************************************************************** *****/
//
//  Y_c[]:               continuous float responses
//                       * ordering: response 0 for cluster 0, response 0 for cluster 1, ..., response 0 for cluster I-1
//                                   response 1 for cluster 0, response 1 for cluster 1, ..., response 1 for cluster I-1
//                                   ...
//                                   response R_c-1 for cluster 0, ..., response R_c-1 for cluster I-1
//
//  Y_d[]:                   discrete integer response
//
//  keepChain_nMCMC_R_cd_dist[]:
//                           0: keepChain = 0/1 indicating whether we should keep chains of "optional" parameters
//                              keepChain[0] = keep sampled values of random effects?
//                           +1: nMCMC[4] = length of MCMC, see NMix_MCMC.h
//                           +4: R_c = number of continuous response variables
//                           +1: R_d = number of discrete response variables
//                               * in the following R = R_c + R_d
//                           +1: dist[R]:                 distribution of each response
//                               * see GLMM::_GLMM_dist in GLMM.h
//
//  I_n[1 + R*I]:            I_n[0] = I = number of groups of correlated observations
//                           I_n[1,...] = n =
//                             numbers of observations for each response and each cluster
//                             * ordering: n[response 0, cluster 0], n[response 0, cluster 1], ..., n[response 0, cluster I-1],
//                                         n[response 1, cluster 0], n[response 1, cluster 1], ..., n[response 1, cluster I-1],
//                                         ...
//                                         n[response R-1, cluster 0], n[response R-1, cluster 1], ..., n[response R-1, cluster I-1], 
//                             * in the following N = sum(n) (N = total number of observations)
//
//  X[]                      covariate matrices for fixed effects (without a column of ones for intercept)
//                           * ordering of covariates: response 0 for observation (0, 0), response 0 for observation (0, 1), 
//                                                     ..., 
//                                                     response 0 for observation (I-1, n[0, I-1]),
//                                                     ...,
//                                                     response R-1 for observation (0, 0), response R-1 for observation (0, 1), 
//                                                     ...,
//                                                     response R-1 for observation (I-1, n[0, I-1])
//  XtX[]                    lower triangles of matrices t(X_s) %*% X_s, where X_s is the design matrix of the fixed effects 
//                           (including possibly column of ones for a fixed intercept) for response s (s=0,...,R-1)
//                           THIS ARGUMENT IS NO MORE PRESENT, XtX matrices are automatically calculated inside the GLMM_MCMC function
//
//  Z[]:                     covariate matrices for random effects (without a column of ones for intercept)
//                           * ordering the same as for X matrix
//
//  SZitZiS[]:               INPUT:  lower triangles of matrices t(Z_s[i]) %*% Z_s[i],
//                                   where Z_s[i] is the design matrix of the random effects
//                                   (including possibly column of ones for a random intercept), sorted in this way:
//                                   t(Z_0[0]) %*% Z_0[0], ..., t(Z_{R-1}[0]) %*% Z_{R-1}[0],
//                                   ...
//                                   t(Z_0[I-1]) %*% Z_0[I-1], ..., t(Z_{R-1}[I-1]) %*% Z_{R-1}[I-1], 
//                           OUTPUT: lower triangles of matrices S_s %*% t(Z_s[i]) %*% Z_s[i] %*% S_s,
//                                   where S is the diagonal matrix with scale_b on a diagonal  
//                           THIS ARGUMENT IS NO MORE PRESENT, SZitZiS matrices are automatically calculated inside the GLMM_MCMC function
//
//  p_fI_q_rI[4*R]:          p_fI_q_rI[0,...,R-1]     = p[R]: 
//                                       numbers of fixed effects covariates for each response corresponding to X matrices
//                                       * these numbers do not count possible intercepts
//                           p_fI_q_rI[R,...,2*R-1]   = fixedIntcpt[R]: 
//                                       0/1 indicating whether there is intercept included in beta vector for each response
//                                       * in the following l_beta = sum(p) + sum(fixedIntcpt)
//                           p_fI_q_rI[2*R,...,3*R-1] = q[R]:
//                                       numbers of random effect covariates for each response corresponding to Z matrices
//                                       * these numbers do not count random intercept
//                           p_fI_q_rI[3*R,...,4*R-1] = randIntcpt[R]:             
//                                       0/1 indicating whether there is random intercept included for each response
//                                       * in the following dim_b = sum(q) + sum(randIntcpt)
//
//  shiftScale_b[2*dim_b]:  shift and scale used for random effects to improve mixing
//
////
//  PARAMETERS OF THE PRIOR DISTRIBUTION
//  ====================================
//
//  priorDouble_eps[]:   double prior (hyper)parameters for distribution of residuals
//                       * priorDouble[+R_c]      = zeta
//                                          * prior degrees of freedom of the Wishart prior on Sigma[j]^{-1} for each response
//                       * priorDouble[+R_c]      = g
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1} for each response
//                       * priorDouble[+R_c]      = h
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1} for each response
//
//  priorInt_b[]:        integer prior (hyper)parameters for distribution of random effects b
//                       * priorInt_b[0]  = type of prior for K (see enum _NMix_type_priorK in NMix.h)
//                                          * 0 -> K is fixed
//                                          * 1 -> uniform on 1, ..., Kmax
//                                          * 2 -> Poisson(lambda) truncated on Kmax
//                       * priorInt_b[+1] = type of prior for (mu[j], Sigma[j]^{-1}) (see enum _NMix_type_priormuSigma in NMix.h)
//                                          * 0 -> natural conjugate
//                                          * 1 -> independent conjugate
//                       * priorInt_b[+1] = Kmax
//                                          (if K is fixed, this is that fixed value of K)
//                       
//
//  priorDouble_b[]:     double prior (hyper)parameters for distribution of random effects b
//                       * priorDouble[0]  = lambda 
//                                     * expectation of the Poisson distribution if the truncated Poisson is used as a prior for K
//                       * priorDouble[+1] = delta 
//                                     * prior "sample size" for the Dirichlet prior on weights
//                       * priorDouble[+dimb*Kmax] = xi[j] (j=0, Kmax-1) 
//                                          * prior means for the mixture means
//                       * priorDouble[+Kmax]   = c[j] (j=0, Kmax-1)
//                                          * prior scale parameters for the mixture means if natural conjugate prior is used
//                       * priorDouble[+LT(dimb)*Kmax]  = D[j]^{-1} (j=0, ..., Kmax-1)                 
//                                          * prior inverse variances for the mixture means if independent conjugate prior is used
//                       * priorDouble[+1]      = zeta
//                                          * prior degrees of freedom of the Wishart prior on Sigma[j]^{-1}
//                       * priorDouble[+dimb]      = g
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1}
//                       * priorDouble[+dimb]      = h
//                                          * parameters of gamma hyperpriors above Sigma[j]^{-1}
//
//  priorDouble_beta[]:  prior (hyper)parameters for distribution of fixed effects beta
//                       * prior_beta[l_beta]  = prior means for beta's
//                       * prior_beta[+l_beta] = prior inverse variances for beta's
//
//
//  PARAMETERS TO TUNE MCMC
//  ================================================================
//
//  tune_scale_beta[R_d]:  scale parameters for each DISCRETE response profile by which we multiply
//                         the proposal covariance matrix when updating the 'fixed' effects
//                         of DISCRETE response profiles
//
//  tune_scale_b[1]:       scale paramater by which we multiply the proposal covariance matrix when updating
//                         the fixed effects
//                         * used only when there are some discrete response profiles
//                         * ignored when there are only continuous responses since then the Gibbs move is used
//                           to update random effects
//
//
//  INITIAL AND LAST SAMPLED VALUES RELATED TO THE DISTRIBUTION OF Y
//  ================================================================
//  sigma_eps[R_c]    INPUT:  initial values for standard deviations of residuals of each continuous response
//                   OUTPUT:  last sampled value
//
//  gammaInv_eps[R_c]
//
//
//  INITIAL AND LAST SAMPLED VALUES RELATED TO THE DISTRIBUTION OF RANDOM EFFECTS b
//  ===============================================================================
//  K_b[1]:                   INPUT:  initial value for number of mixture components in the distribution of random effects
//                           OUTPUT:  last sampled value
//  
//  w_b[Kmax_b]:              INPUT:  initial values for mixture weights in the distribution of random effects
//                           OUTPUT:  last sampled value
//
//  mu_b[Kmax_b*q]:           INPUT:  initial values for mixture means in the distribution of random effects
//                           OUTPUT:  last sampled value
//
//  Q_b[Kmax_b*LT(q)]:        INPUT:  WHATSEVER
//                           OUTPUT:  last sampled values of inverted mixture covariance matrices in the distribution of random effects
//
//  Sigma_b[Kmax_b*LT(q)]:     INPUT: WHATSEVER
//                           OUTPUT:  last sampled values of mixture covariance matrices in the distribution of random effects
//
//  Li_b[Kmax_b*LT(q)]:       INPUT:  initial values for Cholesky decompositions of inverted mixture covariance matrices in the distribution of random effects
//                           OUTPUT:  last sampled value
//
//  gammaInv_b[q]:            INPUT:  initial value of the inverse of the gamma hyperparameter for the distribution of random effects
//                           OUTPUT:  last sampled value of the inverse of the gamma hyperparameter for the distruibution of random effects
//
//  r_b[I]:                   INPUT:  initial allocations for random effects
//                           OUTPUT:  last values of allocations
//
//  r_b_first[I]              INPUT:  whatsever
//                           OUTPUT:  value of r_b corresponding to the first kept (not burn-in, after thinning loop) MCMC iteration
//
//
//  INITIAL AND LAST SAMPLED VALUES RELATED TO REGRESSION
//  =====================================================
//  beta[l_beta]:               INPUT:  initial values of genuine fixed effects
//                             OUTPUT:  last sampled value
//
//  b[dim_b*I]:                 INPUT:  initial values of random effects
//                             OUTPUT:  last sampled value
//
//
//  b_first[dim_b*I]:           INPUT:  whatsever
//                             OUTPUT:  value of b corresponding to the first kept (not burn-in, after thinning loop) MCMC iteration
//
//
//  CHAINS FOR QUANTITIES RELATED TO THE DISTRIBUTION OF RESIDUALS
//  ==============================================================
//  chsigma_eps[]
//
//  chgammaInv_eps[]
//
//  
//  CHAINS FOR QUANTITIES RELATED TO THE DISTRIBUTION OF RANDOM EFFECTS
//  ===================================================================
//  chK_b[]
// 
//  chw_b[]
//
//  chmu_b[]
//
//  chQ_b[]
//
//  chSigma_b[]
//
//  chLi_b[]
//
//  chgammaInv_b[]
//
//  chorder_b[]
//
//  chrank_b[]
//
//  chMeanData_b[]
//
//  chCorrData_b[]
//
//
//  CHAINS FOR REGRESSION RELATED QUANTITIES
//  ========================================
//  chbeta
//
//  chb
//
//  
//  OTHER CHAINS 
//  =============
//  chGLMMLogL[]               chain with the marginal (random effects integrated out) log-likelihood at each iteration of MCMC
//                             - the marginal log-likelihood is obtained using the Laplacian approximation
//
//  chLogL[]                   chain with the conditional (given random effects) log-likelihood at each iteration of MCMC
//
//
//  PERFORMANCE OF MCMC
//  ====================
//  naccept_beta[R_c + R_d]    numbers of accepted beta's (in Metropolis-Hastings step)
//                             for each response profiles
//                             REMARK:  for continuous responses, naccept_beta[s] = number of MCMC iterations
//                                      since Gibbs step is used there
//
//  naccept_b[I]               number of accepted b's (in Metropolis-Hastings step) for each cluster 
//                             REMARK:  for models with continuous response profiles only, naccept_b[i] = number of MCMC iterations
//                                      since Gibbs step is used there
//  
//
//  POSTERIOR MEANS OF SOME QUANTITIES, OTHER POSTERIOR QUANTITIES
//  ==============================================================
//  pm_eta_fixed
//
//  pm_eta_random
//
//  pm_meanY
//
//  pm_stres
//
//  pm_b
//
//  pm_w
//
//  pm_mu_b
//
//  pm_Q_b
//
//  pm_Sigma_b
//
//  pm_Li_b
//
//  pm_indGLMMLogL[I]: posterior means of marginal (random effects integrated out) log-likelihood for each cluster
//
//  pm_indLogL[I]:     posterior means of log-(conditional given b) likelihood for each cluster
//
//  pm_indLogpb[I]:    posterior means of log-density of random effects for each cluster
//
//  sum_Ir_b[I, K_b]:  for each cluster and each mixture component: sum(r_b[i] = k),
//                     i = 0, ..., I-1, j = 0, ..., K_b - 1
//                     from the main part of MCMC (burn-in not included) 
//                     * COMPUTED ONLY WHEN K_b is FIXED
//                     * COMPONENTS ARE (INTERNALLY) RE-LABELED BEFORE sum(r_b[i] = k) is computed
//
//  sum_Pr_b_b[I, K_b]:  for each cluster and each mixture component: sum(P(r_b[i] = k | theta, b, y))
//                       i = 0, ..., I-1, j = 0, ..., K_b - 1
//                       from the main part of MCMC (burn-in not included) 
//                       * COMPUTED ONLY WHEN K_b is FIXED
//                       * COMPONENTS ARE (INTERNALLY) RE-LABELED BEFORE sum(P(r_b[i] = k | theta, b, y)) is computed
//
//  MISCALLANEOUS
//  =============
//  iter[1]:           iteration index
//  
//  err[1]:            error flag
//
void
GLMM_MCMC(double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
          int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
          const int*    keepChain_nMCMC_R_cd_dist,  
          int*          I_n,                                // this is in fact const, not declared as const to be able to use **
          const double* X, 
          double*       Z,                                  // this is in fact const, not declared as const to be able to use **
          const int*    p_fI_q_rI,
          const double* shiftScale_b,
          const double* priorDouble_eps,
          const int*    priorInt_b,           
          const double* priorDouble_b,
          const double* priorDouble_beta, 
          const double* tune_scale_beta,
          const double* tune_scale_b,
          double* sigma_eps,     
          double* gammaInv_eps,
          int*    K_b,              
          double* w_b,             
          double* mu_b,    
          double* Q_b,    
          double* Sigma_b,    
          double* Li_b,
          double* gammaInv_b,    
          int*    r_b,
          int*    r_b_first,
          double* beta,          
          double* b, 
          double* b_first,
          double* chsigma_eps,   
          double* chgammaInv_eps,
          int*    chK_b,            
          double* chw_b,           
          double* chmu_b,  
          double* chQ_b,  
          double* chSigma_b,  
          double* chLi_b,
          double* chgammaInv_b,  
          int*    chorder_b,          
          int*    chrank_b,
          double* chMeanData_b,      
          double* chCorrData_b,
          double* chbeta,        
          double* chb,
          double* chGLMMLogL,
          double* chLogL,
          int*    naccept_beta,
          int*    naccept_b,
          double* pm_eta_fixed,
          double* pm_eta_random,
          double* pm_meanY,
          double* pm_stres,
          double* pm_b,
          double* pm_w_b,         
          double* pm_mu_b,        
          double* pm_Q_b,              
          double* pm_Sigma_b,      
          double* pm_Li_b,
          double* pm_indGLMMLogL,
          double* pm_indLogL,
          double* pm_indLogpb,
          int*    sum_Ir_b,
          double* sum_Pr_b_b,
          int*    iter,
          int*    err);

#ifdef __cplusplus
}
#endif

#endif
