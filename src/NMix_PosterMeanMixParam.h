//
//  PURPOSE:   Function to compute posterior means of mixture parameters
//             while taking into account re-labeling of the MCMC output.
//             Implemented for fixed K
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/02/2010, code taken from the part of NMix_MCMC function
//
//  FUNCTIONS:  
//     * NMix::PosterMeanMixParam  09/02/2010
//
// ======================================================================
//
#ifndef _NMIX_POSTER_MEAN_MIX_PARAM_H_
#define _NMIX_POSTER_MEAN_MIX_PARAM_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::PosterMeanMixParam                                                                  *****/
/***** ***************************************************************************************** *****/
//
// pm_w[K]
//
// pm_mu[p, K]
//
// pm_Q[LT(p), K]
//
// pm_Sigma[LT(p), K]
//
// pm_Li[LT(p), K]
//
// K[1]
//
// chw[K, Mkeep]
//
// chmu[p, K, Mkeep]
//
// chQ[LT(p), K, Mkeep]
//
// chSigma[LT(p), K, Mkeep]
//
// chLi[LT(p), K, Mkeep]
//
// chorder[K, Mkeep]
//
// p[1]
//
// Mkeep[1]
//
void
PosterMeanMixParam(double* pm_w,
                   double* pm_mu,
                   double* pm_Q,
                   double* pm_Sigma,
                   double* pm_Li,
                   const int*    K,
                   const double* chw,
                   const double* chmu,
                   const double* chQ,
                   const double* chSigma,
                   const double* chLi,
                   const int*    chorder,
                   const int*    p,
                   const int*    Mkeep);

}    // end of namespace NMix

#endif
