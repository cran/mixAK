//
//  PURPOSE:   GLMM with a normal mixture in the random effects distribution, 
//             computation of the penalized expected deviance
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111125  created
//
//  FUNCTIONS:  
//     * GLMM_PED  20111206:  first version finished  
//
// ====================================================================================================
//
#ifndef _GLMM_PENALIZED_EXPECTED_DEVIANCE_H_
#define _GLMM_PENALIZED_EXPECTED_DEVIANCE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"
#include "AK_BSTAT.h"

#include "NMix_Utils.h"

#include "GLMM_linear_predictors.h"
#include "GLMM_dY_meanY.h"
#include "GLMM_create_ZS.h"
#include "GLMM_Deviance.h"
#include "GLMM_newData.h"
#include "GLMM_eta_fixed_random2eta_meanY.h"


//namespace GLMM{

#ifdef __cplusplus
extern "C" {
#endif

/***** ************************************************************************** *****/
/***** GLMM_PED                                                                   *****/
/***** ************************************************************************** *****/
//
// PED[5]            INPUT:   whatsever
//                   OUTPUT:  PED[0] = Dbar = observed (expected) deviance
//                            PED[1] = p(opt) = total optimism (see p. 526 of Plummer, 2008)
//                                     obtained without reweighting the sample
//                            PED[2] = PED = PED[0] + PED[1]
//                            PED[3] = p(opt) obtained using the importance sampling
//                            PED[4] = PED[0] + PED[3]
//
// pm_indDevObs[I]   INPUT:   whatsever
//                   OUTPUT:  MCMC estimate of observed deviance for each group of clustered observations,
//                            i.e., pm_indDevObs[i] = average of 
//                            0.5 * (-2 * log(L_obs(theta(chain 1))) -2 * log(L_obs(theta(chain 2))))
//                            over iterations,
//                            where L_obs is the observed data likelihood, i.e., with random effects integrated out
//
// pm_indpopt[I]      INPUT:  whatsever
//                   OUTPUT:  MCMC estimate of individual contributions to optimism obtained
//                            without reweighting the sample
//
// pm_windpopt[I]     INPUT:  whatsever
//                   OUTPUT:  MCMC estimate of individual contributions to optimism obtained
//                            using the importance sampling, see Plummer (2008)
//
// invalid_indDevObs[I]  INPUT:  whatsever
//                      OUTPUT:  number of iterations (in both chains) where log(f(y1))=Inf or log(f(y2))=Inf
//                               -> number between 0 and 2*M
//                               
// invalid_indpopt[I]    INPUT:  whatsever
//                      OUTPUT:  number of iterations where indDevObs = Inf or J(theta1, theta2) = +-Inf
//
// invalid_windpopt[I]   INPUT:  whatsever
//                      OUTPUT:  number of iterations where indDevObs = Inf or J(theta1, theta2) = +-Inf
//                               or weight for importance sampling = Inf
//
// sum_ISweight[I]   INPUT:   whatsever
//                   OUTPUT:  sum of relative weights for importance sampling for each observation
//                            (see p. 530 of Plummer, 2008)
//
// chGLMMLogL1[M]    INPUT:   whatsever
//                   OUTPUT:  (again) calculated values of GLMM log-likelihood based on chain 1 (Deviance = -2 * GLMM LogL)
//
// chGLMMLogL2[M]    INPUT:   whatsever
//                   OUTPUT:  (again) calculated values of GLMM log-likelihood based on chain 1 (Deviance = -2 * GLMM LogL)
//
// chGLMMLogL_repl1_ch1[M]  INPUT:  whatsever
//                          OUTPUT: GLMM log-likelihood based on replicates generated according to parameters from chain 1, 
//                                  evaluated at parameters from chain 1 
//
// chGLMMLogL_repl1_ch2[M]  INPUT:  whatsever
//                          OUTPUT: GLMM log-likelihood based on replicates generated according to parameters from chain 1, 
//                                  evaluated at parameters from chain 2
//
// chGLMMLogL_repl2_ch1[M]  INPUT:  whatsever
//                          OUTPUT: GLMM log-likelihood based on replicates generated according to parameters from chain 2, 
//                                  evaluated at parameters from chain 1 
//
// chGLMMLogL_repl2_ch2[M]  INPUT:  whatsever
//                          OUTPUT: GLMM log-likelihood based on replicates generated according to parameters from chain 2, 
//                                  evaluated at parameters from chain 2
//
// err[1]            INPUT:   whatsever
//                   OUTPUT:  error flag
//
// ===========================================================================================================================
//
// Y_c[]           continuous float responses, see GLMM_MCMC.h
//
// Y_d[]           discrete integer response, see GLMM_MCMC.h
//
// R_cd_dist[]     0:  R_c = number of continuous response variables
//                +1:  R_d = number of discrete response variables
//                           * in the following R = R_c + R_d
//                +1:  dist[R]: distribution of each response
//                     GLMM_eta_fixed_random2eta_meanY.h         * see GLMM::_GLMM_dist in GLMM.h
//
//  I_n[1 + R*I]   I_n[0] = I = number of groups of correlated observations
//                 I_n[1,...] = n = numbers of observations for each response and each cluster, see GLMM_MCMC.h
//                                  * in the following N = sum(n) (N = total number of observations)
//
//  X[]            covariate matrices for fixed effects (without a column of ones for intercept), see GLMM_MCMC.h
//                  
//  Z[]            covariate matrices for random effects (without a column of ones for intercept), see GLMM_MCMC.h
//
//  p_fI_q_rI[4*R]        see GLMM_MCMC.h
//
//  distribution_b[1]     assumed distribution of random effects
//                        * see _Random_effect_dist in GLMM.h 
//
//  shiftScale_b[2*dim_b] see GLMM_MCMC.h
//
//
// ===========================================================================================================================
//
//  chsigma_eps1[R_c * M]           Sampled values of standard deviances of continuous responses, chain 1
//
//  chK_b1[M]                       (Sampled) values of the number of mixture components in the distribution of random effects, chain 1
//                                  = K_b1
//
//  chw_b1[K_b1 * M]                Sampled values of mixture weights, chain 1
//
//  chmu_b1[dim_b * K_b1 * M]       Sampled mixture means, chain 1
//
//  chLi_b1[LT(dim_b) * K_b1 * M]   Cholesky factors of sampled inverse mixture variances, chain 1
//
//  chdf_b1[K_b1 * M]               Sampled values of MVT degrees of freedom (chain 1) if it is assumed that b has a MVT distribution
//
//  chbeta1[ * M]                   Sampled values of regression coefficients, chain 1
//
//  bhat1[dim_b * I]                Estimated values of individual random effects (posterior means, ...), chain 1
//
// ===========================================================================================================================
//
//  chsigma_eps2[R_c * M]           Sampled values of standard deviances of continuous responses, chain 2
//
//  chK_b2[M]                       (Sampled) value of the number of mixture components in the distribution of random effects, chain 2
//                                  = K_b2
//
//  chw_b2[K_b2 * M]                Sampled values of mixture weights, chain 2
//
//  chmu_b2[dim_b* K_b2 * M]        Sampled mixture means, chain 2
//
//  chLi_b2[LT(dim_b) * K_b2 * M]   Cholesky factors of sampled inverse mixture variances, chain 2
//
//  chdf_b2[K_b1 * M]               Sampled values of MVT degrees of freedom (chain 2) if it is assumed that b has a MVT distribution
//
//  chbeta2[ * M]                   Sampled values of regression coefficients, chain 2
//
//  bhat2[dim_b * I]                Estimated values of individual random effects (posterior means, ...), chain 2
//
// ===========================================================================================================================
//
//  M[1]       Lengths of the chains
//
//  Dens_ZERO
//
//  EMin       
//
// ===========================================================================================================================
//
void
GLMM_PED(double*       PED,
         double*       pm_indDevObs,    
         double*       pm_indpopt,    
         double*       pm_windpopt,
         int*          invalid_indDevObs,  
         int*          invalid_indpopt,  
         int*          invalid_windpopt,
         double*       sum_ISweight,    
         //double*       ch_ISweight,
         double*       chGLMMLogL1,
         double*       chGLMMLogL2,
         double*       chGLMMLogL_repl1_ch1,
         double*       chGLMMLogL_repl1_ch2,
         double*       chGLMMLogL_repl2_ch1,
         double*       chGLMMLogL_repl2_ch2,
         int*          err,
         double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
         int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
         const int*    R_cd_dist,  
         int*          I_n,                                // this is in fact const, not declared as const to be able to use **
         const double* X, 
         double*       Z,                                  // this is in fact const, not declared as const to be able to use **
         const int*    p_fI_q_rI,
         const int*    distribution_b,
         const double* shiftScale_b,    
         const double* chsigma_eps1,   
         const int*    chK_b1,            
         const double* chw_b1,           
         const double* chmu_b1,
         const double* chLi_b1,
         const double* chQ_b1,
         const double* chdf_b1,
         const double* chbeta1,        
         const double* bhat1,
         const double* chsigma_eps2,   
         const int*    chK_b2,            
         const double* chw_b2,           
         const double* chmu_b2,
         const double* chLi_b2,
         const double* chQ_b2,
         const double* chdf_b2,
         const double* chbeta2,        
         const double* bhat2,
         const int*    M,
         const double* Dens_ZERO,  
         const double* EMin);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace GLMM ***/

#endif

