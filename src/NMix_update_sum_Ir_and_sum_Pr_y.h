//
//  PURPOSE:   Normal mixture model, 
//             computation of 
//             Pr_y[j, i] = w_j * phi(y_i | mu_j, Sigma_j), i=0,...,n-1, j=0,...,K-1
//             and related cumulative sums
//             cum_Pr_y[j, i] = sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l), i=0,...,n-1, j=0,...,K-1
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/02/2010
//
//  FUNCTIONS:  
//     * update_sum_Ir_and_sum_Pr_y  (2 prototypes) 10/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_UPDATE_SUM_IR_AND_SUM_PR_Y_H_
#define _NMIX_UPDATE_SUM_IR_AND_SUM_PR_Y_H_

#include <R.h>

#include "AK_Basic.h"
#include "AK_Utils.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::update_sum_Ir_and_sum_Pr_y                                                          *****/
/***** PROTOTYPE 1                                                                               *****/
/***** * this prototype also calculates Pr_y from cum_Pr_y                                       *****/
/***** ***************************************************************************************** *****/
//
// sum_Ir[K, n]     INPUT:   whatsever
//                  OUTPUT:  sum_Ir[rank[j], i] = INPUT[rank[j], i] + I[r[i]=j]
//
// sum_Pr_y[K, n]   INPUT:   whatsever
//                  OUTPUT:  sum_Pr_y[rank[j], i] = INPUT[rank[j], i] + Pr_y[j, i]
//
// Pr_y[K, n]       INPUT:   whatsever
//                  OUTPUT:  Pr_y[0, i] = cum_Pr_y[0, i]
//                           Pr_y[j, i] = cum_Pr_y[j, i] - cum_Pr_y[j-1, i], j=1,...,K-1
//
// cum_Pr_y[K, n]   
// 
// r[n]             allocations (before re-labeling)
//
// rank[K]          ranks of mixture components under suitable re-labeling
//
// K[1]             number of mixture components
//
// n[1]             sample size
//
void
update_sum_Ir_and_sum_Pr_y(int*    sum_Ir,
                           double* sum_Pr_y,
                           double* Pr_y,
                           const double* cum_Pr_y,
                           const int* r,
                           const int* rank,
                           const int* K,
                           const int* n);


/***** ***************************************************************************************** *****/
/***** NMix::update_sum_Ir_and_sum_Pr_y                                                          *****/
/***** PROTOTYPE 2                                                                               *****/
/***** * this prototype takes Pr_y on input                                                      *****/
/***** ***************************************************************************************** *****/
//
// sum_Ir[K, n]     INPUT:   whatsever
//                  OUTPUT:  sum_Ir[rank[j], i] = INPUT[rank[j], i] + I[r[i]=j]
//
// sum_Pr_y[K, n]   INPUT:   whatsever
//                  OUTPUT:  sum_Pr_y[rank[j], i] = INPUT[rank[j], i] + Pr_y[j, i]
//
// Pr_y[K, n]       given Pr_y[j, i]
//
// r[n]             allocations (before re-labeling)
//
// rank[K]          ranks of mixture components under suitable re-labeling
//
// K[1]             number of mixture components
//
// n[1]             sample size
//
void
update_sum_Ir_and_sum_Pr_y(int*    sum_Ir,
                           double* sum_Pr_y,
                           const double* Pr_y,
                           const int* r,
                           const int* rank,
                           const int* K,
                           const int* n);

}  // end of namespace NMix

#endif
