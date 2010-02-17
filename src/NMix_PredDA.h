//
//  PURPOSE:   Normal mixture model, 
//             discriminant analysis based on a fitted model
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   15/02/2010
//
//  FUNCTIONS:  
//     * NMix_PredDA  15/02/2010:  
//
// ===================================================================================
//
#ifndef _NORMAL_MIXTURE_PREDICTIVE_DISCRIMINANT_ANALYSIS_H_
#define _NORMAL_MIXTURE_PREDICTIVE_DISCRIMINANT_ANALYSIS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "NMix_Utils.h"
#include "NMix_Pr_y_and_cum_Pr_y.h"
#include "NMix_update_sum_Ir_and_sum_Pr_y.h"

#include "NMix_updateCensObs.h"
#include "NMix_updateAlloc.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDA                                                                               *****/
/***** ***************************************************************************************** *****/
//
//  y0[p, n]                    see NMix_MCMC.h
//
//  y1[p, n]                    see NMix_MCMC.h
//
//  censor[p, n]                see NMix_MCMC.h
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
//  chSigma[LT(p), K, keepMCMC] sampled mixture variances
//
//  chLi[LT(p), K, keepMCMC]    sampled Cholesky decompositions of mixture inverse variances
//
//  chrank[K, keepMCMC]         ranks of the mixture components for each sampled mixture
//                              - C indeces from 0, ..., K-1
//
//  y[p, n]                      INPUT:  initial values for censored observations
//                              OUTPUT:  last sampled values of censored observations
//
//  r[n]                         INPUT:   whatsever
//                              OUTPUT:   last sampled values of component allocation
//                                        - C indeces from 0, ..., K-1
//
//  sum_Ir[n, K]:                INPUT:  whatsever
//                              OUTPUT:  for each observation and each mixture component: sum(r[i] = k),
//                                       i = 0, ..., n-1, j = 0, ..., K - 1
//                                       * components are (internally) re-labeled before sum(r[i] = k) is computed
//
//  hatPr_y[n, K]:               INPUT:  whatsever
//                              OUTPUT:  for each observation and each mixture component: (1/M) * sum(P(r[i] = k | theta, y))
//                                       i = 0, ..., n-1, j = 0, ..., K - 1,
//                                       where M is the number of MCMC iterations
//                                       * components are (internally) re-labeled before sum(P(r[i] = k | theta, y)) is computed
//
//  err[1]                       error flag
//
/***** ***************************************************************************************** *****/
void
NMix_PredDA(const double* y0,
            const double* y1,
            const int*    censor,
            const int*    dimy,
            const int*    keepMCMC,
            const int*    info,
            const int*    K,
            const double* chw,
            const double* chmu,
            const double* chSigma,
            const double* chLi,
            const int*    chrank,
            double* y,
            int*    r,
            int*    sum_Ir,
            double* hatPr_y,
            int*    err);

#ifdef __cplusplus
}
#endif

#endif
