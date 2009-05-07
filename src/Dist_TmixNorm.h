//
//  PURPOSE:   Truncated normal mixture 
//             * random numbers
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/11/2008
//
//  FUNCTIONS:  
//     * Dist::rTmixNorm2  09/11/2008:  
//
// ====================================================================================================
//
#ifndef _DIST_TRUNCATED_MIXTURE_NORM_H_
#define _DIST_TRUNCATED_MIXTURE_NORM_H_

#include <R.h>
#include <Rmath.h>

#include "Dist_Discrete.h"
#include "Dist_TNorm.h"

namespace Dist{

/***** ************************************************************************** *****/
/***** Dist::rTmixNorm1                                                           *****/
/***** ************************************************************************** *****/
//
// x[1]      OUTPUT:  sampled value
//
// K[1]      number of mixture components
//
// cumw[K]   cumulative mixture weights
//
// mu[K]     mixture means
//
// sigma[K]  mixture standard deviations
//
// a[1]      truncation limit 1 (see Dist_TNorm.h)
//
// b[1]      truncation limit 2 (see Dist_TNorm.h)
//
// trunc[1]  type of truncation (see Dist_TNorm.h)
//
void
rTmixNorm1(double* x,
           const int* K,     const double* cumw,  const double* mu,  const double* sigma,
           const double* a,  const double* b,     const int* trunc);

}  /** end of namespace Dist **/

#endif

