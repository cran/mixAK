//
//  PURPOSE:   Implementation of methods declared in GLMM_create_ZiS.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/08/2009
//
// ======================================================================
//
#include "GLMM_create_ZiS.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::create_ZiS                                                                          *****/
/***** ***************************************************************************************** *****/
void
create_ZiS(double* ZiS,     double** ZrespP,
           double** Zresp,  const double* scale_b,  const int* q,  const int* randIntcpt,
           const int* R,    const int* I,           const int* n)
{
  static int s, i, j, k, l;

  static double *ZiSP;

  static const int *nP;
  static const int *q_s, *randIntcpt_s;
  static double *zP;
  static const double *scale_bP, *scale_b_s;  

  /*** Init for some pointers ***/
  for (s = 0; s < *R; s++){
    ZrespP[s] = Zresp[s];
  }

  /*** Loop over longitudinal profiles ***/
  nP   = n;
  ZiSP = ZiS;
  for (i = 0; i < *I; i++){
  
    /*** Loop over observations within a longitudinal profile ***/
    for (j = 0; j < *nP; j++){

      /*** Loop over response types ***/
      q_s          = q;
      randIntcpt_s = randIntcpt;
      scale_b_s    = scale_b;
      for (s = 0; s < *R; s++){

        /*** Loop over rows or resulting Zi_s matrix     ***/
        /*** only up to j                                ***/
        zP = ZrespP[s];
        for (k = 0; k <= j; k++){

          scale_bP    = scale_b_s;

          if (*randIntcpt_s){
            *ZiSP = *scale_bP;
            ZiSP++;
            scale_bP++;            
          }

          for (l = 0; l < *q_s; l++){
            *ZiSP = *scale_bP * *zP;
            ZiSP++;
            zP++;
            scale_bP++;
          }

          if (j == *nP - 1){
            ZrespP[s] = zP;
          }          
        }

        scale_b_s = scale_bP;
        q_s++;
        randIntcpt_s++;
      }        
    }

    nP++;
  }

  return;
}

}
