//
//  PURPOSE:   Implementation of methods declared in Misc_findIndexOfPermutation.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/02/2010
//
// ======================================================================
//
#include "Misc_findIndexOfPermutation.h"

namespace Misc{

/***** ***************************************************************************************** *****/
/***** Misc::findIndexOfPermutation                                                              *****/
/***** ***************************************************************************************** *****/
void
findIndexOfPermutation(int* index,
                       const int* order,
                       const int* order_perm,
                       const int* K,
                       const int* M)
{
  int i, j, row;
  bool notFound;
 
  int *index_now;
  const int *order_now, *orderP, *order_permP;
  
  index_now = index;
  order_now = order;
  for (i = 0; i < *M; i++){      /*** loop over the number of supplied permutations ***/

    notFound = true;
    order_permP = order_perm;
    row = 0;
    while (notFound){

      orderP = order_now;
      j = 0;
      while (j < *K){
        if (*orderP == *order_permP){
          j++;
          orderP++;
          order_permP++;
        }
        else{
          order_permP += (*K - j);
          row++;
          break;
        }
      }

      if (j == *K){
        *index_now = row;
        notFound = false;
      }
    }

    index_now++;
    order_now = orderP;    
  }

  return;
}

}    // end of namespace Misc 
