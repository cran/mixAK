//
//  PURPOSE:   Implementation of methods declared in NMix_reorder_Pr_y.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2010
//
// =============================================================================
//
#include "NMix_reorder_Pr_y.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::reorder_Pr_y                                                                        *****/
/***** ***************************************************************************************** *****/
void
reorder_Pr_y(double *Pr_y,
             double *work,
             const int* order,
             const int* M,
             const int* n,
             const int* K)
{
  int m, i, k;
  const int *orderP;
  
  const int *orderIt = order;
  double *Pr_ySubj = Pr_y;

  for (m = 0; m < *M; m++){
    for (i = 0; i < *n; i++){
      AK_Basic::copyArray(work, Pr_ySubj, *K);
      orderP = orderIt;
      for (k = 0; k < *K; k++){
        *Pr_ySubj = work[*orderP];
        Pr_ySubj++;
        orderP++;
      } 
    }
    orderIt = orderP;
  }

  return;
}

}
