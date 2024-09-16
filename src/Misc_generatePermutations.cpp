//
//  PURPOSE:   Implementation of methods declared in Misc_generatePermutations.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
// ======================================================================
//
#include "Misc_generatePermutations.h"

namespace Misc{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Misc::generatePermutations                                                                *****/
/***** ***************************************************************************************** *****/
void
generatePermutations(int* n_perm,
                     int* order,
                     int* tmp_order,
                     int* rank,
                     const int* K)
{
  int k, j, i, h;
  int *orderP, *tmp_orderP, *tmp_order_now, *rankP;

  if (*K <= 0) Rf_error("Misc::generatePermutations: non-positive K supplied.\n");

  /*** K = 1 ***/
  if (*K == 1){
    *n_perm = 1;
    *order = 0;
    *rank  = 0;
    return;
  }

  /*** K > 1 ***/
  *order  = 0;     /** all permutations of (0)             **/
  *n_perm  = 1;    /** current number of all permutations  **/

  for (k = 1; k < *K; k++){
    
    /*** Copy permutations of (0, ..., k-1), stored in order from the previous loop, to tmp_order ***/
    orderP     = order;
    tmp_orderP = tmp_order;
    for (j = 0; j < *n_perm * k; j++){
      *tmp_orderP = *orderP;
      tmp_orderP++;
      orderP++;
    }

    /*** Merge k to all places of all previous premutations and creates permutations of (0, ..., k) ***/
    orderP = order;
    tmp_order_now = tmp_order;    
    for (j = 0; j < *n_perm; j++){    /*** loop over permutations of (0, ..., k-1)   ***/
      
      /*** new permutation which arises by putting k in front of tmp_order_now ***/
      tmp_orderP = tmp_order_now;
      *orderP = k;
      orderP++;
      for (h = 0; h < k; h++){
        *orderP = *tmp_orderP;
        orderP++;
        tmp_orderP++;
      }

      /*** new permutations which arises by putting k after each element of tmp_order_now  ***/
      for (i = 0; i < k; i++){     /*** loop over elements of tmp_order_now ***/
        tmp_orderP = tmp_order_now;
        for (h = 0; h <= i; h++){
          *orderP = *tmp_orderP;
          orderP++;
          tmp_orderP++;          
        }
        *orderP = k;
        orderP++;
        for (h = i+1; h < k; h++){
          *orderP = *tmp_orderP;
          orderP++;
          tmp_orderP++;          
        }        
      }
   
      tmp_order_now += k;
    }
    

    /*** New number of all permutations ***/
    *n_perm *= (k + 1);
  }

  /*** Calculate 'rank' matrix ("inverse" function to 'order') ***/
  orderP = order;
  rankP  = rank;
  for (j = 0; j < *n_perm; j++){

    for (k = 0; k < *K; k++){
      rankP[*orderP] = k;
      orderP++;
    }  

    rankP += *K;
  }

  return;
}

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Stat ***/

