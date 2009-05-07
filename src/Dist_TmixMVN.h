//
//  PURPOSE:   Truncated multivariate normal mixture 
//             * random numbers
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/11/2008
//
//  FUNCTIONS:  
//     * Dist::TmixMVN2  09/11/2008:  
//
// ====================================================================================================
//
#ifndef _DIST_TRUNCATED_MULTIVARIATE_MIXTURE_NORM_H_
#define _DIST_TRUNCATED_MULTIVARIATE_MIXTURE_NORM_H_

#include <R.h>
#include <Rmath.h>

#include "Dist_Discrete.h"
#include "Dist_TMVN.h"

namespace Dist{

/***** ************************************************************************************ *****/
/***** Dist::rTmixMVN1                                                                      *****/
/***** ************************************************************************************ *****/
//
// x[p]              OUTPUT:  sampled value
//
// K[1]              number of mixture components
//
// cumw[K]           cumulative mixture weights
//
// beta[K * pxp]     beta arguments for Dist::rTMVN1 for each component
// 
// sigmaR2[K * p]    sigmaR2 arguments for Dist::rTMVN1 for each component
//
// a[p]              truncation limit 1 (see Dist_TNorm.h)
//
// b[p]              truncation limit 2 (see Dist_TNorm.h)
//
// trunc[1]          type of truncation (see Dist_TNorm.h)
//
// p[1]              dimension
//
// p_p[1]            p*p
//
void
rTmixMVN1(double* x,  
          const int* K,        const double* cumw,  
          const double* beta,  const double* sigmaR2,  
          const double* a,     const double* b,        const int* trunc,  
          const int* p,        const int* p_p);

}  /** end of namespace Dist **/

#endif

