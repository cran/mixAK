//
//  PURPOSE:   Implementation of methods declared in NMix_Stephens_costMatrix.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/02/2010
//
// =============================================================================
//
#include "NMix_Stephens_costMatrix.h"

namespace NMix{

//#ifdef __cplusplus
//extern "C" {
//#endif

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_costMatrix                                                                 *****/
/***** ***************************************************************************************** *****/
void
Stephens_costMatrix(double* cost,
                    const double* hatPr_y,
                    const double* Pr_y,
                    const int*    n,
                    const int*    K)
{
  static int j, l, i;
  static double *costP;
  static const double *hatPr_yP, *hatPr_y_col_j, *Pr_yP, *Pr_y_col_l;

  costP         = cost;
  Pr_y_col_l    = Pr_y;
  for (l = 0; l < *K; l++){      /*** loop(l) over columns of the cost matrix ***/

    hatPr_y_col_j = hatPr_y;

    for (j = 0; j < *K; j++){      /*** loop(j) over rows of the cost matrix ***/

      *costP   = 0.0;
      hatPr_yP = hatPr_y_col_j;
      Pr_yP    = Pr_y_col_l;

      for (i = 0; i < *n; i++){      /*** loop(i) over subjects ***/

        if (*Pr_yP > AK_Basic::_ZERO){
          if (*hatPr_yP > AK_Basic::_ZERO){ 
            *costP += *Pr_yP * (log(*Pr_yP) - log(*hatPr_yP));
          }
          else{                /*** q_{i,j} = 0, we have to add log(Inf)                ***/
            *costP += 710.0;   /*** in R, last non-Inf log is log(1e309) \approx 710    ***/
          }
        }
        //else{                /*** p_{i,l} = 0, we have to add zero                    ***/
	  //  *costP += 0.0;
        //}

        Pr_yP    += *K;
        hatPr_yP += *K;
      }				     /*** loop(i) end ***/
   
      costP++;
      hatPr_y_col_j++;
    }                              /*** loop(j) end ***/ 

    Pr_y_col_l++;
  }                              /*** loop(l) end ***/

  return;
}

//#ifdef __cplusplus
//}
//#endif

}  // end of namespace NMix
