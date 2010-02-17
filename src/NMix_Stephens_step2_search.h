//
//  PURPOSE:   Normal mixture model, re-labeling algorithm
//
//             * Step 2 of the algorithm given in Stephens (2000, JRSS-B, p. 802)
//
//             * For each iteration m, find such re-labeling (determined by 'order' and 'rank') which minimizes
//               sum_{i=1}^n sum_{j=1}^K Pr_y[m, order[m, j], i] * log(Pr_y[m, order[m, j], i] / hatPr_y[j, i])
//               = sum_{i=1}^n sum_{j=1}^K Pr_y[m, j, i] * log(Pr_y[m, j, i] / hatPr_y[rank[m, j], i])
//             
//             WARNING:  This implementation of step 2 of Stephens' algorithm examines all K! possibilities to find a minimum.
//                       Hence, it is not recommended to call this function with K larger than 5 or 6.
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
//  FUNCTIONS:  
//     * Stephens_step2_search  11/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_STEPHENS_STEP_TWO_SEARCH_H_
#define _NMIX_STEPHENS_STEP_TWO_SEARCH_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step2_search                                                               *****/
/***** ***************************************************************************************** *****/
//
//  nchange[1]         INPUT:  whatsever
//                    OUTPUT:  number of changes in index during this step
//                             nchange = 0 on output indicates that a fixed point has been reached
//
//  index[M]           INPUT:  permutation index vector which corresponds to current solution
//                    OUTPUT:  permutation index vector which corresponds to a new solution
// 
//  order[M, K]        INPUT:  whatsever
//                    OUTPUT:  'order' vectors (for each iteration) which corresponds to optimal solution
// 
//  rank[M, K]         INPUT:  whatsever
//                    OUTPUT:  'rank' vectors (for each iteration) which corresponds to optimal solution
//
//  hatPr_y[n, K]     current values of cluster probabilities for each subject
//                    (corresponding to current re-labeling)
//
//  Pr_y[M, n, K]     cluster probabilities for each iteration and each subject
//                    (before any re-labeling)
//
//  order_perm[K!, K] list of all possible permutations
//
//  M[1]              number of MCMC iterations
//
//  n[1]              number of subjects
//
//  K[1]              number of mixture components
//
//  n_perm[1]         this should be equal to K!
//
// *******************************************************************************************************
void
Stephens_step2_search(int* nchange,
                      int* index,
                      int* order,
                      int* rank,
                      const double* hatPr_y,
                      const double* Pr_y,
                      const int*    order_perm,
                      const int*    M,
                      const int*    n,
                      const int*    K,
                      const int*    n_perm);

}  /*** end of namespace NMix ***/

#endif
