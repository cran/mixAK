//
//  PURPOSE:   Normal mixture model, re-labelling algorithm
//             
//             * computation of the cost matrix for the transportation problem
//               which can be used to perform step 2 of the Stephens' algorithm 
//               (Stephens, 2000, JRSS-B, p. 802)
//
//             WARNING:  Currently, step 2 of Stephens' algorithm is implemented such that all K! possibilities
//                       are examined. Hence, it is not recommended to call this function with K larger than 5 or 6.
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   16/02/2010
//
//  FUNCTIONS:  
//     * Stephens_costMatrix  16/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_STEPHENS_COST_MATRIX_H_
#define _NMIX_STEPHENS_COST_MATRIX_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_costMatrix                                                                 *****/
/***** ***************************************************************************************** *****/
//
//  cost[K, K]           INPUT:  whatsever
//                      OUTPUT:  calculated cost matrix (stored in column major order)
//
//  hatPr_y[n, K]       current values of cluster probabilities (estimated from all MCMC iterations) for each subject
//                      (corresponding to current re-labelling)
//                      * on p. 808 of Stephens (2000), q_{i,j} = hatPr_y[i, j] (hatPr_y is stored in ROW major order when using this notation)
//
//  Pr_y[n, K]          cluster probabilities for a specific MCMC iteration and each subject
//                      (before any re-labelling)
//                      * on p. 808 of Stephens (2000), p_{i,l} = Pr_y[i, l] (Pr_y is stored in ROW major order when using this notation)
//
//  n[1]                number of subjects
//
//  K[1]                number of mixture components
//
// *******************************************************************************************************
void
Stephens_costMatrix(double* cost,
                    const double* hatPr_y,
                    const double* Pr_y,
                    const int*    n,
                    const int*    K);

#ifdef __cplusplus
}
#endif

}    // end of namespace NMix

#endif
