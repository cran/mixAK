//
//  PURPOSE:   Normal mixture model, computation of the predictive conditional densities
//             (all margins given one margin)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   31/05/2009
//
//  FUNCTIONS:  
//     * NMix_PredCondDensMarg  31/05/2009:  
//                        
//
// ====================================================================================================
//
#ifndef _NMIX_PREDICTIVE_CONDITIONAL_DENSITY_MARGINAL_H_
#define _NMIX_PREDICTIVE_CONDITIONAL_DENSITY_MARGINAL_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "Dist_MVN.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredCondDensMarg                                                                     *****/
/***** ***************************************************************************************** *****/
//
// Compute predictive  (univariate) conditional densities for all margins given one chosen margin
//
// IMPLEMENTED ONLY FOR MODELS WITH A FIXED NUMBER OF MIXTURE COMPONENTS
//
// dens[n[icond] + (n[0] + ... + 0 + ... + n[p-1])*n[icond]]
//     OUTPUT:  Computed predictive conditional densities and marginal densities for the margin by which we condition
//              that is: f(ycond[0]), ... f(ycond[n[icond]-1])                                             ..... n[icond] elements
//                       f(y0[0]|ycond=ycond[0]), ..., f(y0[n[0]-1]|ycond=ycond[0])                        ..... n[0] elements
//                       .....
//                       f(y0[0]|ycond=ycond[n[icond]-1]), ..., f(y0[n[0]-1]|ycond=ycond[n[icond]-1])      ..... n[0] elements 
//                       .....
//                       .....
//                       f(yp_1[0]|ycond=ycond[0]), ..., f(yp_1[n[p-1]-1]|ycond=ycond[0])                    ..... n[p-1] elements
//                       .....
//                       f(yp_1[0]|ycond=ycond[n[icond]-1]), ..., f(yp_1[n[p-1]-1]|ycond=ycond[n[icond]-1])  ..... n[p-1] elements 
//
// dwork[2 + LT(p) + n[icond] + (n[0] + ... + 0 + ... + n[p-1])*n[icond]]
//
// err[1]                   Error flag
//
// icond[1] Index of the margin by which we condition
//
// y[n[0] + ... + n[p-1]]   Marginal grids of values to evaluate the conditional densities
//                          The space occuied by the margin by which we condition specifies
//                          the values by which we condition
//                          
// p[1]     Dimension of the distribution
//
// n[p]     Lengths of the marginal grids
//
// chK[1]                Number of mixture components
//
// chw[sum(K)]           Sampled weights
//
// chmu[p*sum(K)]        Sampled means
//
// chLi[LT(p)*sum(K)]    Cholesky factors of sampled inverse variances
//
// M[1]       Lengths of the chains
//
void
NMix_PredCondDensMarg(double* dens,
                      double* dwork,     int* err,
                      const int* icond,  const double* y,  const int* p,       const int* n,  
                      const int* chK,    const double* chw,  const double* chmu,  const double* chLi,
                      const int* M);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif
