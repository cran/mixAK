//
//  PURPOSE:   Normal mixture model, computation of the chains for quantities
//             derived from mixture parameters
//             * marginal means of exp(Y) (useful when we model density of Y = log(T) and T = exp(Y) is of primary interest)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/06/2009
//
//  FUNCTIONS:  
//     * NMix_PredDensMarg  03/12/2007:  
//                          01/02/2008:  Possible scaling added
//
// ====================================================================================================
//
#ifndef _NMIX_CHAINS_DERIVED_H_
#define _NMIX_CHAINS_DERIVED_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix_ChainsDerived                                                                        *****/
/***** ***************************************************************************************** *****/
//
//
// chEexpY[p*M]   
//     OUTPUT:  Chains for mean of exp(Y) (possible shift and scale are taken into account) 
//              stored in this way:
//              EexpY[0,...,p-1] ..... marginal means of exp(Y) at the first iteration
//              EexpY[p,...2*p-1] .... marginal means of exp(Y) at the second iteration
//              ...
//
// dwork[LT(p)]
//
// err[1]                Error flag
//
// p[1]                  Dimension of the distribution
//
// shiftScale[p+p]       shift[p] and scale[p] parameters that should be taken into account
//                       It is assumed that Y = shift + diag(scale) * Z, 
//                                          Z ~ sum_{k=1}^K w_k N(mu_k, Sigma_k)
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
// M[1]                  Lengths of the chains
//
// Krandom[1]            If = 0 then it is assumed that K is constant
//
void
NMix_ChainsDerived(double* chEexpY,
                   double* dwork,    int* err,
                   const int* p,     const double* shiftScale,
                   const int* chK,   const double* chw,         const double* chmu,  const double* chLi,
                   const int* M,     const int* Krandom);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif
