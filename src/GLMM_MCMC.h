//
//  PURPOSE:   Main functions to run MCMC to estimate
//             a GLMM
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/07/2009
//
//  FUNCTIONS:  
//     * GLMM_MCMC  03/07/2009:  Start working on it
//                  03/08/2009:  Version allowing for continuous responses working
//
// =================================================================================
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

#include "GLMM.h"
#include "GLMM_linear_predictors.h"
#include "GLMM_scale_ZitZi.h"
#include "GLMM_updateVars_eps.h"
#include "GLMM_updateHyperVars_eps.h"
#include "GLMM_updateFixEf_gauss.h"
#include "GLMM_updateRanEf_nmix_gauss.h"


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
//  R_c[1]:                  number of continuous response variables
//
//  Y_d[]:                   discrete integer response
//
//  R_d[1]:                  number of discrete response variables
//                           * in the following R = R_c + R_d
//
//  dist[R]:                 distribution of each response
//                           * see GLMM::_GLMM_dist in GLMM.h
//
//  I[1]:                    number of clusters
//
//  n[R*I]:                  numbers of observations for each response and each cluster
//                           * ordering: n[response 0, cluster 0], n[response 0, cluster 1], ..., n[response 0, cluster I-1],
//                                       n[response 1, cluster 0], n[response 1, cluster 1], ..., n[response 1, cluster I-1],
//                                       ...
//                                       n[response R-1, cluster 0], n[response R-1, cluster 1], ..., n[response R-1, cluster I-1], 
//                           * in the following N = sum(n) (N = total number of observations)
//
//  X[]                      covariate matrices for fixed effects (without a column of ones for intercept)
//                           * ordering of covariates: response 0 for observation (0, 0), response 0 for observation (0, 1), ..., response 0 for observation (I-1, n[0, I-1]),
//                                                     ...
//  XtX[]                    lower triangles of matrices t(X_s) %*% X_s, where X_s is the design matrix of the fixed effects 
//                           (including possibly column of ones for a fixed intercept) for response s (s=0,...,R-1)
//                                                     
//  p[R]:                    numbers of fixed effects covariates for each response corresponding to X matrices
//                           * these numbers do not count possible intercepts
//
//  fixedIntcpt[R]:          0/1 indicating whether there is intercept included in beta vector for each response
//                           * in the following l_beta = sum(p) + sum(fixedIntcpt)
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
//                           OUTPUT: lower triangles of matrices S %*% t(Z_s[i]) %*% Z_s[i] %*% S,
//                                   where S is the diagonal matrix with scale_b on a diagonal  
//
//  q[R]:                    numbers of random effect covariates for each response corresponding to Z matrices
//                           * these numbers do not count random intercept
//
//  randIntcpt[R]:           0/1 indicating whether there is random intercept included for each response
//                           * in the following dim_b = sum(q) + sum(randIntcpt)
//
//  shiftScale_b[2*sum(q)]:  shift and scale used for random effects to improve mixing
//
//
//  MCMC PARAMETERS
//  =============== 
//  nMCMC[4]:            length of MCMC, see NMix_MCMC.h
//
//  keepChain[]:         0/1 indicating whether we should keep chains of "optional" parameters
//                       keepChain[0] = keep sampled values of random effects?
//
//
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
//  K_b[1]:                  INPUT:  initial value for number of mixture components in the distribution of random effects
//                          OUTPUT:  last sampled value
//                                   * ignored if distribution is NORMAL
//  
//  w_b[Kmax_b]              INPUT:  initial values for mixture weights in the distribution of random effects
//                          OUTPUT:  last sampled value
//                                   * ignored if distribution is NORMAL
//
//  mu_b[Kmax_b*q]           INPUT:  initial values for mixture means in the distribution of random effects
//                          OUTPUT:  last sampled value
//
//  Q_b[Kmax_b*LT(q)]        INPUT:  WHATSEVER
//                          OUTPUT:  last sampled values of inverted mixture covariance matrices in the distribution of random effects
//
//  Sigma_b[Kmax_b*LT(q)]    INPUT:  WHATSEVER
//                          OUTPUT:  last sampled values of mixture covariance matrices in the distribution of random effects
//
//  Li_b[Kmax_b*LT(q)]       INPUT:  initial values for Cholesky decompositions of inverted mixture covariance matrices in the distribution of random effects
//                          OUTPUT:  last sampled value
//
//  gammaInv_b[q]
//
//  r_b[I]
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
//  POSTERIOR MEANS OF SOME QUANTITIES
//  ==================================
//  pm_eta_fixed
//
//  pm_eta_random
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
//  pm_indLogL[I]:     posterior means of log-(conditional given b) likelihood for each cluster
//
//  pm_indLogpb[I]:    posterior means of log-density of random effects for each cluster
//
//
//  MISCALLANEOUS
//  =============
//  iter[1]:           iteration index
//  
//  err[1]:            error flag
//
void
GLMM_MCMC(double* Y_c,                            // this is in fact const, not declared as const to be able to use **
          const int* R_c,   
          int* Y_d,                               // this is in fact const, not declared as const to be able to use **
          const int* R_d,  
          const int* dist,                 
          const int* I,                  
          int* n,                                 // this is in fact const, not declared as const to be able to use **
          const double* X, 
          const double* XtX,                
          const int* p,                  
          const int* fixedIntcpt,
          double* Z,                              // this is in fact const, not declared as const to be able to use **
          double* SZitZiS,               
          const int* q,                  
          const int* randIntcpt,   
          const double* shiftScale_b,
          const int* nMCMC,
          const int* keepChain,
          const double* priorDouble_eps,
          const int* priorInt_b,           
          const double* priorDouble_b,
          const double* priorDouble_beta, 
          double* sigma_eps,     
          double* gammaInv_eps,
          int* K_b,              
          double* w_b,             
          double* mu_b,    
          double* Q_b,    
          double* Sigma_b,    
          double* Li_b,
          double* gammaInv_b,    
          int* r_b,
          double* beta,          
          double* b, 
          double* chsigma_eps,   
          double* chgammaInv_eps,
          int* chK_b,            
          double* chw_b,           
          double* chmu_b,  
          double* chQ_b,  
          double* chSigma_b,  
          double* chLi_b,
          double* chgammaInv_b,  
          int* chorder_b,          
          int* chrank_b,
          double* chMeanData_b,      
          double* chCorrData_b,
          double* chbeta,        
          double* chb,
          double* pm_eta_fixed,
          double* pm_eta_random,
          double* pm_b,
          double* pm_w_b,         
          double* pm_mu_b,        
          double* pm_Q_b,              
          double* pm_Sigma_b,      
          double* pm_Li_b,
          double* pm_indLogL,
          double* pm_indLogpb,
          int* iter,
          int* err);

#ifdef __cplusplus
}
#endif

#endif
