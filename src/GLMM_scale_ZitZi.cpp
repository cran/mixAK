//
//  PURPOSE:   Implementation of methods declared in GLMM_scale_ZitZi.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/07/2009
//
// ======================================================================
//
#include "GLMM_scale_ZitZi.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::scale_ZitZi                                                                         *****/
/***** ***************************************************************************************** *****/
void
scale_ZitZi(double* SZitZiS,
            const double* scale_b,
            const int* q_ri,        const int* R,  const int* I)
{
  int i, s, j, k;
  double *SZitZiSP = SZitZiS;
  const double *scale_bP  = NULL;
  const double *scale_b2P = NULL;
  const int *q_riP        = NULL;

  for (i = 0; i < *I; i++){
    scale_bP = scale_b;
    q_riP    = q_ri;
    for (s = 0; s < *R; s++){
      for (k = 0; k < *q_riP; k++){    /** loop over columns of Z'Z **/
        scale_b2P = scale_bP;
        for (j = k; j < *q_riP; j++){
          *SZitZiSP *= (*scale_bP * *scale_b2P);
          scale_b2P++;
          SZitZiSP++;
        }
        scale_bP++;
      }    
      q_riP++;
    }
  }
  return;
}

}
