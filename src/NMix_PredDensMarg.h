//
//  PURPOSE:   Normal mixture model, computation of the predictive marginal densities
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/12/2007
//
//  FUNCTIONS:  
//     * NMix_PredDensMarg  03/12/2007:  
//                          01/02/2008:  Possible scaling added
//
// ====================================================================================================
//
#ifndef _NMIX_PREDICTIVE_DENSITY_MARGINAL_H_
#define _NMIX_PREDICTIVE_DENSITY_MARGINAL_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDensMarg                                                                         *****/
/***** ***************************************************************************************** *****/
//
// Compute predictive marginal (univariate) densities for all margins
//
// dens[n[0] + n[1] + ... + n[p-1]]
//     OUTPUT:  Computed predictive marginal densities
//
// densK[(n[0] + n[1] + ... + n[p-1]) * Kmax]
//     OUTPUT:  Computed predictive marginal densities conditioned by K
//              * computed only if Krandom <> 0
//              * if Krandom == 0, NULL can be supplied 
//
// freqK[Kmax]              OUTPUT:  Frequencies of the K values in the sample
//                                   * computed only if Krandom <> 0
//                                   * if Krandom == 0, NULL can be supplied 
//
// propK[Kmax]              OUTPUT:  Posterior distribution of K
//                                   * computed only if Krandom <> 0
//                                   * if Krandom == 0, NULL can be supplied 
//
// dwork[LT(p)]
//
// err[1]                   Error flag
//
// y[n[0] + ... + n[p-1]]   Marginal grids of values to evaluate the density
//
// p[1]     Dimension of the distribution
//
// n[p]     Lengths of the marginal grids
//
// chK[1 or M]           Sampled numbers of mixture components
//                       chK[0] = K if Krandom == 0
//
// chw[sum(K)]           Sampled weights
//
// chmu[p*sum(K)]        Sampled means
//
// chLi[LT(p)*sum(K)]    Cholesky factors of sampled inverse variances
//
// M[1]       Lengths of the chains
//
// Kmax[1]    Maximal K
//
// Krandom[1]  If = 0 then it is assumed that K is constant
//
void
NMix_PredDensMarg(double* dens,     double* densK,   int* freqK,   double* propK,      
                  double* dwork,    int* err,
                  const double* y,  const int* p,       const int* n,  
                  const int* chK,   const double* chw,  const double* chmu,  const double* chLi,
                  const int* M,     const int* Kmax,    const int* Krandom);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif
