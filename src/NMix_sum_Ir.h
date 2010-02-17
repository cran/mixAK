//
//  PURPOSE:   Normal mixture model, 
//             calculate sum_{m=1}^M I[r_i^{(m)} = order(j)] for i=0,...,n-1, j=0,...,K-1
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/02/2010
//
//  FUNCTIONS:  
//     * sum_Ir 14/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_SUM_IR_H_
#define _NMIX_SUM_IR_H_

#include <R.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::sum_Ir                                                                              *****/
/***** ***************************************************************************************** *****///
//  sum_Ir[K, n]      INPUT:   whatsever
//                    OUTPUT:  sum_Ir[rank[j], i] = sum_{m=1}^M I[r^{(m)}[i] = j]
//                             sum_Ir[j, i]       = sum_{m=1}^M I[r^{(m)}[i] = order[j]]
//
//  r[n, M]           allocations (using indeces from the raw MCMC output) 
//
//  rank[K, M]        ranks of mixtire components under suitable re-labeling
//
//  K[1]              number of mixture components
//
//  n[1]              number of subjects
//
//  M[1]              number of MCMC iterations
//
/***** ***************************************************************************************** *****/
void
sum_Ir(int* sum_Ir,
       const int* r,
       const int* rank,
       const int* K,
       const int* n,
       const int* M);

}    // end of namespace NMix

#endif

