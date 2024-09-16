//
//  PURPOSE:   For given K, generate all permutations in the form of the
//             index vector (length K!),
//             order matrix (K! x K),
//             rank matrix (K! x K)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
//  FUNCTIONS:
//     * generatePermutations  11/02/2010
//                    
// =================================================================================================================
//
#ifndef _MISC_GENERATE_PERMUTATIONS_H_
#define _MISC_GENERATE_PERMUTATIONS_H_

#include <R.h>
#include <R_ext/Error.h>

namespace Misc{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Misc::generatePermutations                                                                *****/
/***** ***************************************************************************************** *****/
//
//  n_perm[1]          INPUT:  whatsever
//                    OUTPUT:  number of generated permutations
//
//  order[K!, K]       INPUT:  whatsever
//                    OUTPUT:  permutations of (0,...,K-1)
//
//  tmp_order[K!, K]  working space
// 
//  rank[K!, K]        INPUT:  whatsever
//                    OUTPUT:  'rank' which corresponds to 'order'
//                             order[i, rank[i,j]] = j
//  K[1]           
//
// ****************************************************************************************************
void
generatePermutations(int* n_perm,
                     int* order,
                     int* tmp_order,
                     int* rank,
                     const int* K);

#ifdef __cplusplus
}
#endif

}    /*** end of namespace Misc ***/

#endif
