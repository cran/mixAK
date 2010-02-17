//
//  PURPOSE:   Implementation of methods declared in NMix_Stephens_step1.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/02/2010
//
// =============================================================================
//
#include "NMix_Stephens_step1.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step1                                                                      *****/
/***** ***************************************************************************************** *****/
void
Stephens_step1(double* hatPr_y,
               const double* Pr_y,
               const int*    rank,
               const int*    M,
               const int*    n,
               const int*    K)
{
  static const int *rankP;
  static const double *Pr_yP;

  static double *hatPr_yP;
  static int m, j, i;

  /***** Initialize hatPr_y by zeros *****/
  AK_Basic::fillArray(hatPr_y, 0.0, *n * *K);

  /***** Loop over MCMC iterations to calculate sum_{m=0}^{M-1} Pr_y[m, order(j), i] *****/
  rankP = rank;
  Pr_yP = Pr_y;
  for (m = 0; m < *M; m++){
    
    hatPr_yP = hatPr_y;
    for (i = 0; i < *n; i++){

      for (j = 0; j < *K; j++){
        hatPr_yP[rankP[j]] += *Pr_yP;
        Pr_yP++;
      }
      hatPr_yP += *K;     /** jump to the next subject        **/
    
    }
    rankP += *K;          /** jump to the next MCMC iteration **/
  }

  /***** Calculate averages over MCMC iterations - final estimates *****/
  hatPr_yP = hatPr_y;
  for (i = 0; i < *n * *K; i++){
    *hatPr_yP /= *M;
    hatPr_yP++;
  }

  return;
}

}    /*** end of namespace NMix ***/

