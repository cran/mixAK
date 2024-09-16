//
//  PURPOSE:   Wishart distribution
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/11/2007
//             * some functions partially taken from Mvtdist3.{h,cpp} of the glmmAK package 
//
//  FUNCTIONS:
//     * rWishartEye      12/11/2007:   Sample from the Wishart distribution W(nu, Eye)
//                                      * (almost) completely taken from rwishartEye3 function in Mvtdist3.{h,cpp}[glmmAK]
//
//     * rWishart         12/11/2007:   Sample from the Wishart distribution W(nu, S)
//                                      * partially taken from rwishart3 function in Mvtdist3.{h,cpp}[glmmAK]
//
//     * rWishart_diagS   28/01/2008:   Sample from the Wishart distribution W(nu, S) when S is diagonal
//
//     * l_Wishart_const  24/01/2008:   Logarithm of the factor in the Wishart density which depends only on degrees of freedom
//
//     * ldWishart0       15/01/2008:   Log-density of the Wishart distribution (more basic version)
//
//     * ldWishart        15/01/2008:   Log-density of the Wishart distribution
//
//     * ldWishart_diagS  16/01/2008:   Log-density of the Wishart distribution with diagonal scale matrix
//
//     * rWishart_R       12/11/2007:   R wrapper for rWishart
//
//     * ldWishart_R      15/01/2008:   R wrapper for ldWishart0 and ldWishart
//
// ==============================================================================================================================
//
#ifndef _DIST_WISHART_H_
#define _DIST_WISHART_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>       
#include <R_ext/Error.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rWishartEye                                                                         *****/
/***** ***************************************************************************************** *****/
/*                                                                                   */
/* rwishartEye: Sample from the Wishart distribution W(nu, Eye),                     */
/*              where Eye stands for an identity matrix                              */
/*                                                                                   */
/* Bartlett's decomposition as described on p. 99 of                                 */
/* Ripley (1987). Stochastic simulation. New York: Wiley is used.                    */
/*                                                                                   */
/* See also p. 41 of my red notes.                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*    W[LT(dim)]       OUTPUT:  sampled matrix                                       */
/*                                                                                   */
/*    dwork[LT(dim)]   working array of length LT(dim)                               */
/*                                                                                   */
/*    nu[1]:           degrees of freedom of the Wishart distribution                */
/*                                                                                   */
/*    dim[1]:          dimension of the Wishart distribution                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
void
rWishartEye(double* W, double* dwork, const double* nu, const int* dim);


/***** ***************************************************************************************** *****/
/***** Dist::rWishart                                                                            *****/
/***** ***************************************************************************************** *****/
/*                                                                                   */
/* rWishart: Sample from the Wishart distribution W(nu, S)                           */
/*                        (parametrization as in Gelman, 2004),                      */
/*                        i.e., expectation is nu*S                                  */ 
/*                                                                                   */
/* Algorithm from Ripley (1987), p. 99 is used.                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*            Let S^{-1} = Li * t(Li)                                                */
/*                     S = t(Li^{-1}) * Li^{-1}                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*    W[LT(dim)]      OUTPUT:  sampled matrix                                        */
/*                                                                                   */
/*    dwork[2*dim^2]  working array of length 2*dim^2                                */
/*                                                                                   */
/*    nu[1]           degrees of freedom of the Wishart distribution                 */
/*                                                                                   */
/*    Li[LT(dim)]    lower triangle of the matrix dim x dim with Li,                 */
/*                   where S^{-1} = Li*t(Li)                                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
void
rWishart(double* W,  double* dwork,  const double* nu,  const double* Li,  const int* dim);


/***** ***************************************************************************************** *****/
/***** Dist::rWishart_diagS                                                                      *****/
/***** ***************************************************************************************** *****/
/*                                                                                   */
/* rWishart: Sample from the Wishart distribution W(nu, S)                           */
/*           with diagonal scale matrix S                                            */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */ 
/*    W[LT(dim)]       OUTPUT:  sampled matrix                                       */
/*                                                                                   */
/*    dwork[LT(dim)]   working array of length LT(dim)                               */
/*                                                                                   */
/*    nu[1]:           degrees of freedom of the Wishart distribution                */
/*                                                                                   */
/*    d_invS[dim]:     diagonal of S^{-1}                                            */
/*                                                                                   */
/*    dim[1]:          dimension of the Wishart distribution                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
void
rWishart_diagS(double* W,  double* dwork,  const double* nu,  const double* d_invS,  const int* dim);


/***** ***************************************************************************************** *****/
/***** Dist::l_Wishart_const                                                                     *****/
/***** ***************************************************************************************** *****/
//
// Computes logarithm of the factor in the Wishart density which only depends on degrees of freedom
//
// log_const[1]         Logarithm of the factor in the density which only depends on degrees of freedom
//
// nu[1]                Degrees of freedom of the Wishart distribution
//
// dim[1]               Dimension of the Wishart distribution
//
void
l_Wishart_const(double* log_const,  const double* nu,   const int* dim);


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart0                                                                          *****/
/***** Dist::ldWishart                                                                           *****/
/***** ***************************************************************************************** *****/
//
// log_dens[1]          Computed log-density of the Wishart distribution
//
// W[LT(dim)]           Lower triangle of the matrix for which we want to evaluate a Wishart density
//
// W_L[LT(dim)]         Cholesky decomposition of W
//
// log_sqrt_detW[1]     log|W|^{1/2} = 0.5*log|W|
//
// log_const[1]         Logarithm of the factor in the density which only depends on degrees of freedom
//
// nu[1]                Degrees of freedom of the Wishart distribution
//
// invS[LT(dim)]        Lower triangle of the inverted scale matrix
//
// invS_L[LT(dim)]      Cholesky decomposition of invS
//
// log_sqrt_detinvS[1]  log|S|^{-1/2} = -0.5*log|S| = 0.5*log|S|^{-1}
//
// dim[1]               Dimension of the Wishart distribution
//
void 
ldWishart0(double* log_dens,  double* log_sqrt_detW,  double* log_const,     double* log_sqrt_detinvS,  
           const double* W,   const double* W_L,      
           const double* nu,  const double* invS,     const double* invS_L,  const int* dim);

void
ldWishart(double* log_dens,
          const double* W,          const double* log_sqrt_detW, 
          const double* log_const,  const double* nu,   
          const double *invS,       const double* log_sqrt_detinvS,  const int* dim);


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart_diagS                                                                     *****/
/***** ***************************************************************************************** *****/
//
// Compute log-density of the Wishart distribution when the scale matrix S is diagonal.
//
// log_dens[1]          Computed log-density of the Wishart distribution
//
// W[LT(dim)]           Lower triangle of the matrix for which we want to evaluate a Wishart density
//
// log_sqrt_detW[1]     log|W|^{1/2} = 0.5*log|W|
//
// log_const[1]         Logarithm of the factor in the density which only depends on degrees of freedom
//
// nu[1]                Degrees of freedom of the Wishart distribution
//
// invS_diag[dim]       Diagonal elements of the inverted scale matrix
//
// log_sqrt_detinvS[1]  log|S|^{-1/2} = -0.5*log|S| = 0.5*log|S|^{-1}
//
// dim[1]               Dimension of the Wishart distribution
//
void
ldWishart_diagS(double* log_dens,
                const double* W,          const double* log_sqrt_detW, 
                const double* log_const,  const double* nu,   
                const double *invS_diag,  const double* log_sqrt_detinvS,  const int* dim);


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rWishart_R                                                                          *****/
/***** ***************************************************************************************** *****/
//
// R wrapper to rWishart
//
// x[LT(dim)*npoints]:   OUTPUT: sampled values (lower triangles only)
//
// dwork[2*dim^2]:       working array
//
// nu[1]:                Wishart degrees of freedom
//
// invS[LT(dim)]:        S^{-1} = inverse scale matrix of the Wishart distribution
//                       * expectation is nu * S
//
// dim:                  dimension of the Wishart distribution
//
// npoints:              number of sampled points
//
void
rWishart_R(double* W,  double* dwork,  int* err,  const double* nu,  double* invS,  const int* dim,  const int* npoints);


/***** ***************************************************************************************** *****/
/***** Dist::ldWishart_R                                                                         *****/
/***** ***************************************************************************************** *****/
//
// log_dens[npoints]:
//
// W_L[npoints*LT(dim)]:    OUTPUT:   Cholesky decompositions of matrices W
//
// log_sqrt_detW[npoints]:  
//
// log_const[1]:            
//
// invS_L[LT(dim)]:         OUTPUT:   Cholesky decomposition of invS 
//                      
// log_sqrt_detinvS[1]:     
//
// err[1]:
//
// W[npoints*LT(dim)]:      Input matrices W (their lower triangles only)
// 
// nu[1]:
//
// invS[LT(dim)]:
// 
// dim[1]:
//
// npoints[1]:
//
void 
ldWishart_R(double* log_dens,   double* W_L,       double* log_sqrt_detW,  
            double* log_const,  double* invS_L,    double* log_sqrt_detinvS,  int* err,
            const double* W,    const double* nu,  const double* invS,        
            const int* dim,     const int* npoints);

#ifdef __cplusplus
}
#endif

}   /*** end of namespace Dist ***/

#endif

