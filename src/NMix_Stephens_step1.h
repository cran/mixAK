//
//  PURPOSE:   Normal mixture model, re-labeling algorithm
//
//             * Step 1 of the algorithm given in Stephens (2000, JRSS-B, p. 802)
//
//             * Computation of 
//             hatPr_y[j, i] = (1 / M) * sum_{m=1}^M Pr_y[m, order[m, j], i],
//             where  
//             Pr_y[m, j, i] = w_j^{(m)} * phi(y_i | mu_j^{(m)}, Sigma_j^{(m)}), i=0,...,n-1, j=0,...,K-1, m=0,...,M-1
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
//  FUNCTIONS:  
//     * Stephens_step1  11/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_STEPHENS_STEP_ONE_H_
#define _NMIX_STEPHENS_STEP_ONE_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step1                                                                      *****/
/***** ***************************************************************************************** *****/
//
//  hatPr_y[n, K]      INPUT: whatsever
//                    OUTPUT: calculated estimated cluster probabilities for each subject
//                            (using supplied re-labeling indicated by 'rank' argument)
//
//  Pr_y[M, n, K]     cluster probabilities for each iteration and each subject
//                    (before any re-labeling)
//
//  rank[M, K]        ranks of components under specific re-labeling
//
//  M[1]              number of MCMC iterations
//
//  n[1]              number of subjects
//
//  K[1]              number of mixture components
//
// *******************************************************************************************************
void
Stephens_step1(double* hatPr_y,
               const double* Pr_y,
               const int*    rank,
               const int*    M,
               const int*    n,
               const int*    K);

}  /*** end of namespace NMix ***/

#endif
