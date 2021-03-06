//
//  PURPOSE:   Normal mixture model, update of the mixture weights
//             and related quantities
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007             
//
//  FUNCTIONS:  
//     * updateWeights  07/11/2007:  Update of the weights
//                      30/03/2015:  Factor covariate on mixture weight allowed
//
// ======================================================================
//
#ifndef _NMIX_UPDATE_WEIGHTS_H_
#define _NMIX_UPDATE_WEIGHTS_H_

#include <R.h>

#include "AK_Basic.h"

#include "Dist_Dirichlet.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateWeights                                                                       *****/
/***** ***************************************************************************************** *****/
//
// Prior:             w       ~ Dirichlet(delta, ..., delta)
// Full conditional:  w | ... ~ Dirichlet(delta+mixN[0], ..., delta+mixN[K-1])
//
// w[K, nxw]            INPUT:   whatsever
//                     OUTPUT:  updated weights
//
// logw[K, nxw]         INPUT:   whatsever
//                     OUTPUT:  updated log-weights
//
// dwork[K]        working array
//
// mixN[K]         numbers of observations belonging to each component
//                 * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// K[1]            current number of components
//
// delta[1]        prior hyperparameter
//
// mixNxw[K, nxw]  numbers of observations belonging into each component
//                 if there is a factor covariate on mixture weights
//                 - added on 20150330
// 
// nxw[1]          number of levels of a factor covariate on mixture weights
//                 - added on 20150330
//
void
updateWeights(double* w,           double* logw,  double* dwork,
              const int* mixN,     const int* K,  
              const double* delta,
              const int* mixNxw,   const int* nxw);

}   /** end of namespace NMix **/

#endif
