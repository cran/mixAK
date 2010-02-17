//
//  PURPOSE:   Normal mixture model, re-labeling algorithm
//
//             * Step 2 of the algorithm given in Stephens (2000, JRSS-B, p. 802)
//
//             * For each iteration m, find such re-labeling (determined by 'order' and 'rank') which minimizes
//               sum_{i=1}^n sum_{j=1}^K Pr_y[m, order[m, j], i] * log(Pr_y[m, order[m, j], i] / hatPr_y[j, i])
//               = sum_{i=1}^n sum_{j=1}^K Pr_y[m, j, i] * log(Pr_y[m, j, i] / hatPr_y[rank[m, j], i])
//          
//             * This implementation uses the corresponding transportation problem to minimize the loss function
//               and should be reasonable fast even for higher values of K
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   16/02/2010
//
//  FUNCTIONS:  
//     * Stephens_step2_transport  16/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_STEPHENS_STEP_TWO_TRANSPORT_H_
#define _NMIX_STEPHENS_STEP_TWO_TRANSPORT_H_

#include <R.h>

#include "AK_Basic.h"

//#include "lpSolve.h"

#include "NMix_Stephens_costMatrix.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step2_transport                                                            *****/
/***** ***************************************************************************************** *****/
//
//  nchange[1]         INPUT:  whatsever
//                    OUTPUT:  number of changes in ordering during this step
//                             nchange = 0 on output indicates that a fixed point has been reached
//
//  order[M, K]        INPUT:  whatsever
//                    OUTPUT:  'order' vectors (for each iteration) which corresponds to optimal solution
// 
//  rank[M, K]         INPUT:  whatsever
//                    OUTPUT:  'rank' vectors (for each iteration) which corresponds to optimal solution
//
//  lp_costs[1 + (K, K)]     INPUT:  whatsever
//                          OUTPUT:  working space to store cost matrix for the transportation problem
//                                   (the first element is equal to 0 and then the costs matrix follows)
//
//  lp_solution[K, K]  INPUT:  whatsever
//                    OUTPUT:  working space to store optimal solution of the transportation problem
//
//
//  The following lp_* variables are not altered by this function and are assumed to be equal to the following on the input:
//  ------------------------------------------------------------------------------------------------------------------------
//
//  lp_r_signs[K]     vector (3, 3, ..., 3) of length K (input for lp_transbig)
//
//  lp_r_rhs[K]       vector (1, 1, ..., 1) of length K (input for lp_transbig)
//
//  lp_c_signs[K]     vector (3, 3, ..., 3) of length K (input for lp_transbig)
//
//  lp_c_rhs[K]       vector (1, 1, ..., 1) of length K (input for lp_transbig)
//
//  lp_integers[K]    vector (1, 2, ..., K) (input for lp_transbig)
//
//  hatPr_y[n, K]     current values of cluster probabilities for each subject
//                    (corresponding to current re-labeling)
//
//  Pr_y[M, n, K]     cluster probabilities for each iteration and each subject
//                    (before any re-labeling)
//
//  M[1]              number of MCMC iterations
//
//  n[1]              number of subjects
//
//  K[1]              number of mixture components
//
/***** ***************************************************************************************** *****/
void
Stephens_step2_transport(int*    nchange,
                         int*    order,
                         int*    rank,
                         double* lp_costs,
                         double* lp_solution,
                         int*    lp_r_signs,
                         double* lp_r_rhs,
                         int*    lp_c_signs,
                         double* lp_c_rhs,
                         int*    lp_integers,
                         const double* hatPr_y,
                         const double* Pr_y,
                         const int*    M,
                         const int*    n,
                         const int*    K);

}    // end of namespace NMix

#endif
