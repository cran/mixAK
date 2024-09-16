//
//  PURPOSE:   (Multivariate) GLMM.
//             Calculate eta   = eta_fixed + eta_random
//                       meanY = invers. link(eta) 
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111206  created
//
//  FUNCTIONS:  
//     * GLMM::eta_fixed_random2eta_meanY  20111206: working version
//
// ====================================================================================================
//
#ifndef _GLMM_ETA_FIXED_RANDOM2ETA_MEANY_H_
#define _GLMM_ETA_FIXED_RANDOM2ETA_MEANY_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"

#include "GLMM.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::eta_fixed_random2eta_meanY                                                          *****/
/***** ***************************************************************************************** *****/
//
//  dist[R]:                      types of response
//
//  n[I, R]:                      number of observations for each cluster and each response
//
//  R[1]:                         number of responses
//
//  I[1]:                         number of clusters
//
void
eta_fixed_random2eta_meanY(double* eta,
                           double* meanY,
                           const double* eta_fixed,
                           const double* eta_random,
                           const int* dist,
                           const int* n,
                           const int* R,
                           const int* I);

}  /*** end of namespace GLMM ***/

#endif


