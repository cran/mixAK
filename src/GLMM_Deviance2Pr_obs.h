//
//  PURPOSE:   Use some output from GLMM::Deviance and calculate P(r[i] = k | theta, y)
//             (random effects not in the condition)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/11/2010
//
//  FUNCTIONS:  
//     * GLMM_Deviance2Pr_obs   29/11/2010
//
// ======================================================================
//
#ifndef _GLMM_DEVIANCE_TO_PR_OBS_H_
#define _GLMM_DEVIANCE_TO_PR_OBS_H_

#include <R.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ********************************************************************** *****/
/***** GLMM:Deviance2Pr_obs                                                   *****/
/***** ********************************************************************** *****/
//
//  Pr_obs[K, I]     INPUT: whatsever
//                  OUTPUT: calculated values of P(u[i] = k | theta, y)
//
//  marg_L_i[I]      values of the individual contributions to the likelihood,
//                   typically output from GLMM::Deviance
//                   which satisfies marg_L_i[i] = sum_k (w[k] * marg_L_ik[k, i]) for all i
// 
//  marg_L_ik[K, I]  values of the individual and component specific contributions to the likelihood,
//                   typically output from GLMM::Deviance
// 
//  w[K]             mixture weights
//
//  I[1]             number of subjects
//
//  K[1]             number of mixture components
//
void
Deviance2Pr_obs(double* Pr_obs,
                const double* marg_L_i,
                const double* marg_L_ik,
                const double* w,
                const int*    I,               
                const int*    K);

}    // end of namespace GLMM

#endif
