//
//  PURPOSE:   Main functions to re-label the MCMC output obtained using the NMix_MCMC function
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/02/2010
//
//  FUNCTIONS:  
//     * NMix_NMixRelabel  10/02/2010:  
//                         29/11/2010:  argument Pr_b_b added and all P(u_i=k | b, theta, y) over all iterations are stored 
//                         31/03/2015:  factor covariate on mixture weights allowed
//
// ======================================================================
//
#ifndef _NORMAL_MIXTURE_NMIX_RELABEL_H_
#define _NORMAL_MIXTURE_NMIX_RELABEL_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

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

#include "NMix_updateCensObs.h"
#include "NMix_updateAlloc.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_NMixRelabel                                                                          *****/
/***** ***************************************************************************************** *****/
//
//  type[1]                     type of the re-labeling algorithm (see enum _NMix_relabel_algorithm in NMix.h)
//
//  iparam[]                    integer parameters for the re-labeling algorithm
//                              type = 1 = NMix::MEAN
//                              =====================
//                                iparam[0] = margin of the mixture means used to determine ordering
//
//                              type = 2 = NMix::WEIGHT
//                              =======================
//                                iparam = NULL
//
//                              type = 3 = NMix::STEPHENS
//                              =========================
//                                iparam[0] = type of the initial re-labeling algorithm
//                                            0 = NMix::IDENTITY
//                                            1 = NMix::MEAN
//                                            2 = NMix::WEIGHT
//                                iparam[1] = margin of the mixture means used to determine initial ordering
//                                            (if initial ordering = NMix::MEAN)
//                                iparam[2] = maximal number of iterations of the re-labeling algorithm
//                                
//  y0[p, n]                    see NMix_MCMC.h
//
//  y1[p, n]                    see NMix_MCMC.h
//
//  censor[p, n]                see NMix_MCMC.h
//
//  nxw_xw[1 + n]   information on a factor covariate on mixture weights (ADDED ON 20150327)
//                  nxw_xw[0]   = nxw = number of covariate levels
//                  nxw_xw[1:n] = covariate values (numbered 0, 1, ..., nxw - 1)
//
//  dimy[2]                     see NMix_MCMC.h
//
//  keepMCMC[1]                 length of the chains
//
//  info[1]                     frequency to report progress of computation
//
//  K[1]                        number of mixture components
//
//  chw[K, keepMCMC]            sampled mixture weights
//
//  chmu[p, K, keepMCMC]        sampled mixture means
//
//  chQ[LT(p), K, keepMCMC]     sampled mixture inverse variances
//
//  chSigma[LT(p), K, keepMCMC] sampled mixture variances
//
//  chLi[LT(p), K, keepMCMC]    sampled Cholesky decompositions of mixture inverse variances
//
//  chorder[K, keepMCMC]         INPUT:  whatsever
//                               OUTPUT: ordering of the mixture components for each sampled mixture
//                                       - C indeces from 0, ..., K-1
//
//  chrank[K, keepMCMC]          INPUT:  whatsever
//                               OUTPUT: ranks of the mixture components for each sampled mixture
//                                       - C indeces from 0, ..., K-1
//
//  y[p, n]                      INPUT:   initial values for censored observations
//                               OUTPUT:  last sampled values of censored observations
//
//  r[n]                         INPUT:   whatsever
//                               OUTPUT:  last sampled values of component allocation
//                                        - C indeces from 0, ..., K-1
//
//  pm_w[K]                      INPUT:  whatsever
//                               OUTPUT: posterior mean of mixture weights
//                                       * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_mu[p, K]                  INPUT:  whatsever
//                               OUTPUT: posterior mean of mixture means
//                                       * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Q[LT(p), K]               INPUT:  whatsever
//                               OUTPUT: posterior mean of mixture inverse variances
//                                       * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Sigma[LT(p), K]           INPUT:  whatsever
//                               OUTPUT: posterior mean of mixture variances
//                                       * before the posterior mean is computed, components are (internally) re-labeled
//
//  pm_Li[LT(p), K]              INPUT:  whatsever
//                               OUTPUT: posterior mean of Cholesky decompositions of mixture inverse variances
//                                       * before the posterior mean is computed, components are (internally) re-labeled
//
//  sum_Ir[K, n]:                INPUT:  whatsever
//                               OUTPUT: for each observation and each mixture component: sum(r[i] = k),
//                                       i = 0, ..., n-1, j = 0, ..., K - 1
//                                       * components are (internally) re-labeled before sum(r[i] = k) is computed
//
//  hatPr_y[K, n]:                INPUT: whatsever
//                               OUTPUT: for each observation and each mixture component: (1/M) * sum(P(r[i] = k | theta, y))
//                                       i = 0, ..., n-1, j = 0, ..., K - 1,
//                                       where M is the number of MCMC iterations
//                                       * components are (internally) re-labeled before sum(P(r[i] = k | theta, y)) is computed
//
//  Pr_y[K, n, keepMCMC]:         INPUT:  whatsever
//                               OUTPUT:  posterior sample of P(r[i] = k | theta, y)
//                                        * columns are re-shuffled to correspond to final re-labelling
//                                                                                                 
//
//  iter_relabel[1]              INPUT:   whatsever
//                               OUTPUT:  unaltered for simple re-labeling algorithms
//                                        For Stephens' algorithm, it gives the number of re-labeling iterations
//
//  nchange[]                    INPUT:   whatsever
//                               OUTPUT:  unaltered for simple re-labeling algorithms
//                                        For Stephens' algorithm, it gives the number of labelling changes at each
//                                        re-labeling iteration.
//                                        * nchange must be of length iparam[2]
//                                        * at convergence nchange[iter_relabel-1] should be equal to 0
//
//  err[1]                       error flag
//
/***** ***************************************************************************************** *****/
void
NMix_NMixRelabel(const int*    type,
                 const int*    iparam,
                 const double* y0,
                 const double* y1,
                 const int*    censor,
                 const int*    nxw_xw,
                 const int*    dimy,
                 const int*    keepMCMC,
                 const int*    info,
                 const int*    K,
                 const double* chw,
                 const double* chmu,
                 const double* chQ,
                 const double* chSigma,
                 const double* chLi,                  
                 int*    chorder,
                 int*    chrank,
                 double* y,
                 int*    r,
                 double* pm_w,
                 double* pm_mu,
                 double* pm_Q,
                 double* pm_Sigma,
                 double* pm_Li,
                 int*    sum_Ir,
                 double* hatPr_y,
                 double* Pr_y,
                 int*    iter_relabel,
                 int*    nchange,
                 int*    err);

#ifdef __cplusplus
}
#endif

#endif
