//
//  PURPOSE:   Implementation of methods declared in NMix_Stephens_step2_search.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
// =============================================================================
//
#include "NMix_Stephens_step2_search.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step2_search                                                               *****/
/***** ***************************************************************************************** *****/
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
                      const int*    n_perm)
{
  static const double *hatPr_yP, *Pr_y_now, *Pr_yP;
  static const int *order_perm_now, *order_permP;

  static int *indexP, *orderP, *rankP;
  static int m, ip, i, j;
  static int index_minLoss;  
  static double minLoss, Loss, pij, pij_q;

  *nchange = 0;
  indexP   = index;
  orderP   = order;
  rankP    = rank;
  Pr_y_now = Pr_y;
  for (m = 0; m < *M; m++){    /*** loop(m) over number of MCMC iterations ***/

    minLoss       = R_PosInf;
    index_minLoss = -1;

    order_perm_now = order_perm;
    for (ip = 0; ip < *n_perm; ip++){   /*** loop(ip) over all possible permutations ***/

      Loss     = 0.0;
      hatPr_yP = hatPr_y;
      Pr_yP    = Pr_y_now;      
      for (i = 0; i < *n; i++){           /*** loop(i) over subjects ***/
        
        order_permP = order_perm_now;
        for (j = 0; j < *K; j++){           /*** loop(j) over mixture components ***/

          pij = Pr_yP[*order_permP];
          if (pij > AK_Basic::_ZERO){
            if (*hatPr_yP > AK_Basic::_ZERO){ 
              Loss += pij * (log(pij) - log(*hatPr_yP));
            }
            else{              /*** q_{i,j} = 0, we have to add log(Inf)              ***/
	      Loss += 710.0;   /*** in R, last non-Inf log is log(1e309) \approx 710  ***/
            }
          }
          //else{              /*** p_{i,j} = 0, we have to add zero                  ***/
	    //  Loss += 0.0;
          //}

          hatPr_yP++;
          order_permP++;
        }                                   /*** end of loop(j) over mixture components ***/

        Pr_yP += *K;
      }                                   /*** end of loop(i) over subjects ***/

      if (Loss < minLoss){                /*** check whether this might be a minimum ***/
        minLoss       = Loss;
        index_minLoss = ip; 
      }

      order_perm_now += *K;
    }                                   /*** end of loop(ip) over all possible permutations ***/

    if (index_minLoss != *indexP){      /*** there is a change in the optimal permutation ***/
      *nchange += 1;
      *indexP = index_minLoss;

      order_permP = order_perm + *K * index_minLoss;
      for (j = 0; j < *K; j++){
        *orderP = *order_permP;
        rankP[*orderP] = j; 

        orderP++;
        order_permP++;
      }
      rankP += *K;
      indexP++;
    }
    else{                              /*** there is no change in the optimal permutation ***/
      indexP++;
      orderP += *K;
      rankP += *K;      
    }

    Pr_y_now = Pr_yP;
  }                            /*** end of loop(m) over number of MCMC iterations ***/

  return;
}

}    /*** end of namespace NMix ***/
