//
//  PURPOSE:   Multivariate normal mixture
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
//  FUNCTIONS:
//     * dmixMVN    06/11/2008:  Multivariate normal mixture distribution, density
//
//     * rmixMVN    07/11/2008:  Multivariate normal mixture distribution, random number + density
//
//     * dmixMVN_R  07/11/2008:  R wrapper to dmixMVN  
//
//     * rmixMVN_R  07/11/2008:  R wrapper to rmixMVN
//                    
// ==============================================================================================================================
//
#ifndef _DIST_MIXTURE_MULTIVARIATE_NORMAL_H_
#define _DIST_MIXTURE_MULTIVARIATE_NORMAL_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"
#include "Dist_Discrete.h"

namespace Dist{

/***** ******************************************************************************** *****/
/***** Dist::dmixMVN                                                                    *****/
/***** ******************************************************************************** *****/
//
//  dens[1]         computed density
//
//  work[nx]        working array
//
//  x[nx]           value where the density is to be evaluated
//
//  K[1]            number of mixture components
//
//  w_dets[K]       mixture w[k] * |Li[k]| * (2*pi)^(-nx/2)
//                  - |Li[k]| = product of diagonal elements of Li[k]
//
//  mu[K*nx]        mixture means
//
//  Li[K*LT(nx)]    Cholesky decompositions of mixture precision matrices
//                  (lower triangles only)
//
//  nx[1]           dimension of the mixture
// 
void
dmixMVN(double* dens,      double* work,
        const double* x,  
        const int* K,      const double* w_dets,  
        const double* mu,  const double* Li,        const int* nx);


/***** ******************************************************************************** *****/
/***** Dist::rmixMVN                                                                    *****/
/***** ******************************************************************************** *****/
//
//  x[nx]           sampled value
//
//  dens[1]         computed density
//
//  work[nx]        working array
//
//  K[1]            number of mixture components
//
//  w_dets[K]       mixture w[k] * |Li[k]| * (2*pi)^(-nx/2)
//                  - |Li[k]| = product of diagonal elements of Li[k]
//
//  cumw[K]         cumulative mixture weights
//
//  mu[K*nx]        mixture means
//
//  Li[K*LT(nx)]    Cholesky decompositions of mixture precision matrices
//                  (lower triangles only)
//
//  nx[1]           dimension of the mixture
// 
void
rmixMVN(double* x,         double* dens,          double* work,
        const int* K,      const double* w_dets,  const double* cumw,  
        const double* mu,  const double* Li,      const int* nx);


/***** ******************************************************************************** *****/
/***** Dist::dmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
//
//  dens[npoints]   computed densities
//
//  w_dets[K]       INPUT:   mixture weights
//                  OUTPUT:  w[k] * |Li[k]| * (2*pi)^(-nx/2) 
//
//  Li[K*LT(nx)]    INPUT:   mixture precision matrices
//                  OUTPUT:  Cholesky decompositions of mixture precision matrices
//
//  work[nx]        working array
//
//  err[1]          error flag
// 
//  x[npoints*nx]   values where the density is to be evaluated
//
//  K[1]            number of mixture components
//
//  mu[K*nx]        mixture means
//
//  nx[1]           dimension of the mixture
//
//  npoints[1]      number of points where the density is to be evaluated
// 
#ifdef __cplusplus
extern "C" {
#endif

void
dmixMVN_R(double* dens,      double* w_dets,   double* Li, 
          double* work,      int* err,
          const double* x,   const int* K,     const double* mu,  
          const int* nx,     const int* npoints);

#ifdef __cplusplus
}
#endif


/***** ******************************************************************************** *****/
/***** Dist::rmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
//
//  x[npoints*nx]   sampled values
//
//  dens[npoints]   computed densities
//
//  w_dets[K]       INPUT:   mixture weights
//                  OUTPUT:  w[k] * |Li[k]| * (2*pi)^(-nx/2) 
//
//  cumw[K]         INPUT:   arbitrary
//                  OUTPUT:  cummulative weights
//
//  Li[K*LT(nx)]    INPUT:   mixture precision matrices
//                  OUTPUT:  Cholesky decompositions of mixture precision matrices
//
//  work[nx]        working array
//
//  err[1]          error flag
// 
//  K[1]            number of mixture components
//
//  mu[K*nx]        mixture means
//
//  nx[1]           dimension of the mixture
//
//  npoints[1]      number of points where the density is to be evaluated
// 
#ifdef __cplusplus
extern "C" {
#endif

void
rmixMVN_R(double* x,         double* dens,     double* w_dets,   double* cumw,  double* Li, 
          double* work,      int* err,
          const int* K,      const double* mu,  
          const int* nx,     const int* npoints);

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Dist ***/

#endif
