//
//  PURPOSE:   Normal mixture model, update of the hyperparameters of the variances
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2007
//
//  FUNCTIONS:  
//     * updateHyperVars  26/11/2007:  Update of the hyperparameters of the variances
//
// ======================================================================
//
#ifndef _NMIX_UPDATE_HYPER_VARS_H_
#define _NMIX_UPDATE_HYPER_VARS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateHyperVars                                                                     *****/
/***** ***************************************************************************************** *****/
//
// gammaInv[p]          INPUT:  whatsever
//                      OUTPUT: updated gamma^{-1}'s
//
// XiInv[LT(p)]         INPUT:  whatsever
//                      OUTPUT: matrix with gamma^{-1}'s on a diagonal
//                              = inverse of the prior scale matrix of the Wishart distribution
//                                * a priori E(Sigma_j^{-1}) = zeta * Xi
//
// log_sqrt_detXiInv[1]   INPUT:  whatsever
//                        OUTPUT: log|XiInv|^{1/2} = 0.5*sum[i]log(gamma[i]^{-1})
//
// dwork[p]             working array
//
// Q[LT(p), K]          inverse mixture variances
//
// K[1]                 current number of components
//
// p[1]                 dimension of the normal distribution
//
// zeta[1]              prior Wishart degrees of freedom
// 
// g[p]                 prior gamma shapes
//
// h[p]                 prior gamma rates
//
void
updateHyperVars(double* gammaInv,    double* XiInv,    double* log_sqrt_detXiInv,  double* dwork,
                const double* Q,     const int* K,     const int* p,
                const double* zeta,  const double* g,  const double* h);

}    /*** end of namespace NMix ***/

#endif

