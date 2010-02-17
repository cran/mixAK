//
//  PURPOSE:   This complements function declared in Misc_generatePermutations.h.
//             It finds the index of the row in the matrix of all permutations which corresponds to supplied permutation.
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/02/2010
//
//  FUNCTIONS:
//     * findIndexOfPermutation 12/02/2010
//                    
// =================================================================================================================
//
#ifndef _MISC_FIND_INDEX_OF_PERMUTATION_H_
#define _MISC_FIND_INDEX_OF_PERMUTATION_H_

#include <R.h>

namespace Misc{

/***** ***************************************************************************************** *****/
/***** Misc::findIndexOfPermutation                                                              *****/
/***** ***************************************************************************************** *****/
//
// index[M]        INPUT:  whatsever
//                OUTPUT:  found permutation indeces for each row of 'order'
//
// order[M, K]         permutations for which we want to find indeces
//
// order_perm[K!, K]   list of all permutations of (0, ..., K-1)
//
// K[1]
//
// M[1]
//
// n_perm[1]           this should be equal to K! (total number of permutations)
//
/***** ***************************************************************************************** *****/
void
findIndexOfPermutation(int* index,
                       const int* order,
                       const int* order_perm,
                       const int* K,
                       const int* M);

}    /*** end of namespace Misc ***/

#endif


