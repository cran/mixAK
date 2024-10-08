//
//  PURPOSE:   Normal mixture model, computation of the predictive pairwise joint densities
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/11/2007
//
//  FUNCTIONS:  
//     * NMix_PredDensJoint2  03/12/2007:  
//
// ====================================================================================================
//
#ifndef _NMIX_PREDICTIVE_DENSITY_JOINT_TWO_H_
#define _NMIX_PREDICTIVE_DENSITY_JOINT_TWO_H_

#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"
#include "Dist_MVN.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_PredDensJoint2                                                                       *****/
/***** ***************************************************************************************** *****/
//
// Compute predictive bivariate densities for all pairs (there are choose(p, 2) = p*(p-1)/2 pairs)
//
// dens[n[0]*n[1] + n[0]*n[2] + ... + n[0]*n[p-1] + ... + n[p-2]*n[p-1]]
//     OUTPUT:  Computed predictive bivariate densities
//
// densK[(n[0]*n[1] + n[0]*n[2] + ... + n[0]*n[p-1] + ... + n[p-2]*n[p-1]) * Kmax]
//     OUTPUT:  Computed predictive bivariate densities conditioned by K
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
// dwork[2+2+2+LT(p)+2+LT(2)]
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
NMix_PredDensJoint2(double* dens,     double* densK,     int* freqK,     double* propK,      
                    double* dwork,    int* err,
                    const double* y,  const int* p,       const int* n,  
                    const int* chK,   const double* chw,  const double* chmu,  const double* chLi,
                    const int* M,     const int* Kmax,    const int* Krandom);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif
