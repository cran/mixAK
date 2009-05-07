//
//  PURPOSE:   Normal mixture model, update of the censored observations.
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   27/11/2007
//
//  FUNCTIONS:  
//     * updateCensObs  27/11/2007:  Update of the censored observations under the normal mixture model
//
// ====================================================================================================
//
#ifndef _NMIX_UPDATE_CENSORED_OBSERVATIONS_H_
#define _NMIX_UPDATE_CENSORED_OBSERVATIONS_H_

#include <R.h>

#include "Stat_BLA.h"
#include "Dist_TNorm.h"
#include "Dist_TMVN.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateCensObs                                                                       *****/
/***** ***************************************************************************************** *****/
//
// It is assumed that there are some censored observations
// (otherwise this function call is wasting of time)
//
// y[p, n]    INPUT:   current values of the censored observations
//            OUTPUT:  updated censored observations
//
// beta[p, p, K]   INPUT:   whatsever
//                 OUTPUT:  if p > 1: regression coefficients for each component (see Stat::BLA)
//
// sigmaR2[p, K]   INPUT:   whatsever
//                 OUTPUT:  if p > 1: residual variances for each component (see Stat::BLA)
//
// dwork[K or LT(p-1)]    when p = 1: working space to store mixture standard deviations
//                        when p > 1: working space of length LT(p-1) for Stat::BLA
//
// err[1]
//
// y0[p, n]
//
// y1[p, n]
//
// censor[p, n]
//
// r[n]
//
// mu[p, K]
//
// Sigma[LT(p), K]
// 
//
// K[1]
//
// p[1]
//
// n[1]
//
void
updateCensObs(double* y,         double* beta,      double* sigmaR2,      double* dwork,  int* err,          
              const double* y0,  const double* y1,  const int* censor,  
              const int* r,      const double* mu,  const double* Sigma,  
              const int* K,      const int* p,      const int* n);

}  /*** end of namespace NMix ***/

#endif
