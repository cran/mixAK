//
//  PURPOSE:   Univariate normal mixture
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
//  FUNCTIONS:
//     * dmixNorm    06/11/2008:     Univariate normal mixture distribution, density
//
//     * rmixNorm    06/11/2008:     Univariate normal mixture distribution, random number + density
//
//     * dmixNorm_R  07/11/2008:     R wrapper to dmixNorm 
//
//     * rmixNorm_R  07/11/2008:     R wrapper to rmixNorm
//                      
// ==============================================================================================================================
//
#ifndef _DIST_MIXTURE_NORM_H_
#define _DIST_MIXTURE_NORM_H_

#include <R.h>
#include <Rmath.h>

#include "Dist_Discrete.h"

namespace Dist{

/***** ********************************************************************************* *****/
/***** Dist::dmixNorm                                                                    *****/
/***** ********************************************************************************* *****/
//
//  dens[1]     computed density
//
//  x[1]        value where the density is to be evaluated
//
//  K[1]        number of mixture components
//
//  w[K]        mixture weights
//
//  mu[K]       mixture means
//
//  sigma[K]    mixture standard deviations
// 
void
dmixNorm(double* dens,
         const double* x,  
         const int* K,  const double* w,  const double* mu,  const double* sigma);


/***** ********************************************************************************* *****/
/***** Dist::rmixNorm                                                                    *****/
/***** ********************************************************************************* *****/
//
//  x[1]        sampled value
//
//  dens[1]     mixture density evaluated in the sampled value
//
//  K[1]        number of mixture components
//
//  w[K]        mixture weights
//
//  cumw[K]     cumulative mixture weights
//
//  mu[K]       mixture means
//
//  sigma[K]    mixture standard deviations
// 
void
rmixNorm(double* x,   double* dens,
         const int* K,  const double* w,  const double* cumw,  const double* mu,  const double* sigma);


/***** ******************************************************************************** *****/
/***** Dist::dmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
//
//  dens[npoints]     INPUT:   arbitrary
//                    OUTPUT:  computed densities
//
//  x[npoinst]        values where the density is to be evaluated
//
//  K[1]              number of mixture components
//
//  w[K]              mixture weights
//
//  mu[K]             mixture means
//
//  sigma[K]          mixture standatd deviations
//
//  npoints[1]        number of points where the density is to be evaluated
// 
#ifdef __cplusplus
extern "C" {
#endif

void
dmixNorm_R(double* dens,
           const double* x,   const int* K,     const double* w,  const double* mu,  const double* sigma,
           const int* npoints);

#ifdef __cplusplus
}
#endif


/***** ******************************************************************************** *****/
/***** Dist::rmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
//
//  x[npoints]        INPUT:   arbitrary
//                    OUTPUT:  sampled points
//
//  dens[npoints]     INPUT:   arbitrary
//                    OUTPUT:  computed densities
//
//  cumw[K]           INPUT:   arbitrary
//                    OUTPUT:  cummulative weights
//
//  K[1]              number of mixture components
//
//  w[K]              mixture weights
//
//  mu[K]             mixture means
//
//  sigma[K]          mixture standatd deviations
//
//  npoints[1]        number of points to be sampled
// 
#ifdef __cplusplus
extern "C" {
#endif

void
rmixNorm_R(double* x,          double* dens,     double* cumw,
           const int* K,       const double* w,  const double* mu,  const double* sigma,  
           const int* npoints);

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Dist ***/

#endif
