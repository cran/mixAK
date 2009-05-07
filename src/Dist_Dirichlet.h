//
//  PURPOSE:   Random number generation from a Dirichlet distribution D(alpha_0,...,alpha_{K-1})
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007
//
//  FUNCTIONS:
//     * rDirichlet  07/11/2007:   Sample from a Dirichlet distribution D(alpha_0,...,alpha_{K-1})
//                                 
//     * rDirichlet_R  07/11/2007:    R wrapper for rDirichlet which allows to generate several vectors points in a loop
//                                    * GetRNGstate() and PutRNGstate() is used inside!!!
//                     
// ==============================================================================================================================
//
#ifndef _DIST_DIRICHLET_H_
#define _DIST_DIRICHLET_H_

#include <R.h>
#include <Rmath.h>

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rDirichlet                                                                          *****/
/***** ***************************************************************************************** *****/
//
// Sample from a Dirichlet distribution D(alpha_0,...,alpha_{K-1}) through sampling from Gamma distribution
// as described in Devroye (1986), p. 564 or in Gelman, Carlin, Stern, Rubin (2004), p. 582.
//
// sampledw[K]      OUTPUT:  sampled vector of probabilities
//
// alpha[K]         parameters of the Dirichlet distribution
//                  (`prior sample sizes')
// 
// K[1]             dimension of the Dirichlet distribution
//
// * No input checks concerning alpha parameters are performed!
//
void
rDirichlet(double* sampledw,  const double* alpha,  const int* K);


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rDirichlet_R                                                                        *****/
/***** ***************************************************************************************** *****/
//
// R wrappers to rDirichlet.
//
// sampledw[K, npoints]   sampled vectors
//
// alpha[K]               the same as in rDirichlet
//
// K[1]                   the same as in rDirichlet
//
// npoints[1]             number of points to sample
//
void
rDirichlet_R(double* sampledw,  const double* alpha,  const int* K,  const int* npoints);

#ifdef __cplusplus
}
#endif

}    /*** end of namespace Dist ***/

#endif
