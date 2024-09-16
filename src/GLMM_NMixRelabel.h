//
//  PURPOSE:   Main functions to re-label the MCMC output obtained using the GLMM_MCMC function
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/02/2010
//
//  FUNCTIONS:  
//     * GLMM_NMixRelabel  26/02/2010:  
//                         27/11/2010:  argument Pr_b_b added and all P(u_i=k | b, theta, y) over all iterations are stored 
//
// ======================================================================
//
#ifndef _GLMM_NMIX_RELABEL_H_
#define _GLMM_NMIX_RELABEL_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Error.h>
#include <R_ext/RS.h>

#include "AK_Basic.h"
#include "AK_BSTAT.h"

#include "Misc_generatePermutations.h"
#include "Misc_findIndexOfPermutation.h"

#include "NMix.h"
#include "NMix_Utils.h"
#include "NMix_orderComp.h"
#include "NMix_PosterMeanMixParam.h"
#include "NMix_Pr_y_and_cum_Pr_y.h"
#include "NMix_update_sum_Ir_and_sum_Pr_y.h"
#include "NMix_sum_Ir.h"
#include "NMix_Stephens_step1.h"
#include "NMix_Stephens_step2_search.h"
#include "NMix_Stephens_step2_transport.h"
#include "NMix_reorder_Pr_y.h"

#include "NMix_updateAlloc.h"
#include "GLMM_updateRanEf.h"
#include "GLMM_updateRanEf_QR.h"

#include "GLMM_linear_predictors.h"
#include "GLMM_create_SZitZiS.h"
#include "GLMM_dY_meanY.h"

#include "GLMM_create_ZS.h"
#include "GLMM_Deviance.h"
#include "GLMM_Deviance2Pr_obs.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_NMixRelabel                                                                          *****/
/***** ***************************************************************************************** *****/
//
//  type[1]                               type of the re-labeling algorithm (see enum _NMix_relabel_algorithm in NMix.h)
//
//  iparam[]                              integer parameters for the re-labeling algorithm
//                                        * see NMix_NMixRelabel.h for details
//
//  nonSilent[1]                          if != 0 then some information is printed during computation
//
//  Y_c[]                                 continuous response 
//                                        * see GLMM_MCMC.h for details
//
//  Y_d[]                                 discrete response 
//                                        * see GLMM_MCMC.h for details
//
//  R_cd[2]                               number of continuous and discrete responses
//                                        * see GLMM_MCMC.h for details
//                                        * in the following, R = R_cd[0] + R_cd[1]
//
//  dist[R]                               types of response
//                                        * see GLMM_MCMC.h for details
//  
//  I[1]                                  number of groups of correlated observations
//                                        * see GLMM_MCMC.h for details
//
//  n[R*I]                                numbers of observations for each response and each cluster
//                                        * see GLMM_MCMC.h for details
//
//  X[]                                   design matrix for fixed effects
//                                        * see GLMM_MCMC.h for details
//
//  Z[]                                   design matrix for random effects
//                                        * see GLMM_MCMC.h for details
//
//  p_fI_q_rI[4]                          dimensions of the model
//                                        * see GLMM_MCMC.h for details
//                                        * in the following, dim_b = dimension of random effects
//
//  shiftScale_b[2*dim_b]                 shift and scale vectors
//                                        * see GLMM_MCMC.h for details
//
//  keepMCMC[1]                           length of the chains
//
//  info[1]                               frequency to report progress of computation
//
//  tune_scale_b[1]                       tuning parameter for update of random effects
//                                        * see GLMM_MCMC.h for details 
//
//  chsigma_eps[R_c, keepMCMC]            sampled values of standard deviations of each continuous response
//
//  distribution_b[1]                     assumed distribution of random effects
//                                        * see _Random_effect_dist in GLMM.h 
//
//  K_b[1]                                number of mixture components
//
//  chw_b[K_b, keepMCMC]                  sampled mixture weights
//
//  chmu_b[dim_b, K_b, keepMCMC]          sampled mixture means
//
//  chQ_b[LT(dim_b), K_b, keepMCMC]       sampled mixture inverse variances
//
//  chSigma_b[LT(dim_b), K_b, keepMCMC]   sampled mixture variances
//
//  chLi_b[LT(dim_b), K_b, keepMCMC]      sampled Cholesky decompositions of mixture inverse variances
//
//  chdf_b[K_b, keepMCMC]                 sampled values of degrees of freedom if random effects are assumed to follow
//                                        the MVT distribution
//
//  chbeta[#fixed effects, keepMCMC]      sampled values of fixed effects
//
//  chorder_b[K_b, keepMCMC]              INPUT:  whatsever
//                                        OUTPUT: ordering of the mixture components for each sampled mixture
//                                                - C indeces from 0, ..., K_b-1
//
//  chrank_b[K_b, keepMCMC]               INPUT:  whatsever
//                                        OUTPUT: ranks of the mixture components for each sampled mixture
//                                                - C indeces from 0, ..., K_b-1
//
//  b[dim_b, I]                           INPUT:   initial values for random effects (on the original scale)
//                                        OUTPUT:  last sampled values of random effects observations
//
//  r_b[I]                                INPUT:   whatsever
//                                        OUTPUT:  last sampled values of component allocation
//                                                 - C indeces from 0, ..., K_b-1
//
//  naccept_b[I]                          INPUT:  whatsever
//                                        OUTPUT: number of accepted b's (in Metropolis-Hastings step) for each cluster 
//                                                * see GLMM_MCMC.h for details
//
//  pm_w_b[K_b]                           INPUT:  whatsever
//                                        OUTPUT: posterior mean of mixture weights
//                                                * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_mu_b[dim_b, K_b]                   INPUT:  whatsever
//                                        OUTPUT: posterior mean of mixture means
//                                                * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Q_b[LT(dim_b), K_b]                INPUT:  whatsever
//                                        OUTPUT: posterior mean of mixture inverse variances
//                                                * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Sigma_b[LT(dim_b), K_b]            INPUT:  whatsever
//                                        OUTPUT: posterior mean of mixture variances
//                                                * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Li_b[LT(dim_b), K_b]               INPUT:  whatsever
//                                        OUTPUT: posterior mean of Cholesky decompositions of mixture inverse variances
//                                                * before the posterior mean is computed, components are (internally) re-labeled
//
//  sum_Ir_b[K_b, I]:                     INPUT:  whatsever
//                                        OUTPUT: for each group of correlated observations and each mixture component: sum(r_b[i] = k),
//                                                i = 0, ..., I-1, j = 0, ..., K_b - 1
//                                                * components are (internally) re-labeled before sum(r_b[i] = k) is computed
//
//  hatPr_b_b[K_b, I]:                    INPUT:  whatsever
//                                        OUTPUT: for each observation and each mixture component: (1/M) * sum(P(r_b[i] = k | theta, b, y))
//                                                i = 0, ..., I-1, j = 0, ..., K_b - 1,
//                                                where M is the number of MCMC iterations
//                                                * components are (internally) re-labeled before sum(P(r_b[i] = k | theta, b, y)) is computed
//
//  Pr_b_b[K_b, I, keepMCMC]:              INPUT:  whatsever
//                                        OUTPUT:  posterior sample of P(r_b[i] = k | theta, b, y)
//                                                 * columns are re-shuffled to correspond to final re-labelling
//                                                                                                 
//  hatPr_obs[K_b, I]:                    INPUT:  whatsever
//                                        OUTPUT: for each observation and each mixture component: (1/M) * sum(P(r_b[i] = k | theta, y)),
//                                                i.e., random effects are (numerically) integrated out
//                                                i = 0, ..., I-1, j = 0, ..., K_b - 1,
//                                                where M is the number of MCMC iterations
//                                                * components are (internally) re-labeled before sum(P(r_b[i] = k | theta, y)) is computed
//
//  Pr_obs[K_b, I, keepMCMC]:              INPUT:  whatsever
//                                        OUTPUT:  posterior sample of P(r_b[i] = k | theta, y),
//                                                 i.e., random effects are (numerically) integrated out
//                                                 * columns are re-shuffled to correspond to final re-labelling
//                                                                                                 
//  iter_relabel[1]                       INPUT:   whatsever
//                                        OUTPUT:  unaltered for simple re-labeling algorithms
//                                                 For Stephens' algorithm, it gives the number of re-labeling iterations
//
//  nchange[]                             INPUT:   whatsever
//                                        OUTPUT:  unaltered for simple re-labeling algorithms
//                                                 For Stephens' algorithm, it gives the number of labelling changes at each
//                                                 re-labeling iteration.
//                                                 * nchange must be of length iparam[2]
//                                                 * at convergence nchange[iter_relabel-1] should be equal to 0
//
//  err[1]                                error flag
//
/***** ***************************************************************************************** *****/
void
GLMM_NMixRelabel(const int*    type,
                 const int*    iparam,
                 const int*    nonSilent,
                 double*       Y_c,                                // this is in fact const, not declared as const to be able to use **
                 int*          Y_d,                                // this is in fact const, not declared as const to be able to use **
                 const int*    R_cd,  
                 const int*    dist,                 
                 const int*    I,                  
                 int*          n,                                  // this is in fact const, not declared as const to be able to use **
                 const double* X, 
                 double*       Z,                                  // this is in fact const, not declared as const to be able to use **
                 const int*    p_fI_q_rI,
                 const double* shiftScale_b,
                 const int*    keepMCMC,
                 const int*    info,
                 const double* tune_scale_b,
                 const double* chsigma_eps,
                 const int*    distribution_b,
                 const int*    K_b,
                 const double* chw_b,
                 const double* chmu_b,
                 const double* chQ_b,
                 const double* chSigma_b,
                 const double* chLi_b,
                 const double* chdf_b,
                 const double* chbeta,                  
                 int*    chorder_b,
                 int*    chrank_b,
                 double* b,
                 int*    r_b,
                 int*    naccept_b,
                 double* pm_w_b,
                 double* pm_mu_b,
                 double* pm_Q_b,
                 double* pm_Sigma_b,
                 double* pm_Li_b,
                 int*    sum_Ir_b,
                 double* hatPr_b_b,
                 double* Pr_b,
                 double* hatPr_obs,
                 double* Pr_obs,
                 int*    iter_relabel,
                 int*    nchange,
                 int*    err);

#ifdef __cplusplus
}
#endif

#endif
