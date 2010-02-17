//
//  PURPOSE:   Implementation of methods declared in NMix_sum_Ir.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/02/2010
//
// =============================================================================
//
#include "NMix_sum_Ir.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::sum_Ir                                                                              *****/
/***** ***************************************************************************************** *****/
void
sum_Ir(int* sum_Ir,
       const int* r,
       const int* rank,
       const int* K,
       const int* n,
       const int* M)
{
  static int l, j, m;
  static int *sum_IrP;
  static const int *rP, *rankP;

  rP    = r;
  rankP = rank;
  for (m = 0; m < *M; m++){      /*** loop over MCMC iterations ***/
     
    sum_IrP   = sum_Ir;
    for (l = 0; l < *n; l++){    /*** loop over subjects ***/
      sum_IrP[rankP[*rP]]++;
      rP++;
      sum_IrP += *K;
    }

    rankP += *K;
  }

  return;
}

}    // end of namespace NMix
