//
//  PURPOSE:   Implementation of methods declared in GLMM_dY_meanY.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/02/2010 as GLMM_dY_d_mean_Y_d.cpp
//             12/04/2010 as this file
//
// ======================================================================
//
#include "GLMM_dY_meanY.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM:dY_mean_Y                                                                            *****/
/***** ***************************************************************************************** *****/
void
dY_meanY(double* dY,
         double* sum_dY_i,
         double* sum_dY,
         double* meanY,
         int*    err,
         const double* Y_c,
         const int*    Y_d,
         const double* eta,
         const int*    dist,
         const int*    n,
         const int*    I,
         const int*    R_c,
         const int*    R_d)
{
  int s, i, j;

  const double *Y_cP  = Y_c;
  const int    *Y_dP  = Y_d;
  const double *etaP  = eta;
  const int    *distP = dist;
  const int    *nP    = n;

  double *dYP    = dY;
  double *meanYP = meanY;

  double *sum_dY_iP;

  AK_Basic::fillArray(sum_dY_i, 0.0, *I);

  for (s = 0; s < (*R_c + *R_d); s++){
    switch (*distP){
    case GLMM::GAUSS_IDENTITY:
      for (i = 0; i < *I; i++){
        for (j = 0; j < *nP; j++){
          *dYP    = 0;
          *meanYP = *etaP;

          Y_cP++;
          dYP++;
          meanYP++;
          etaP++;
        }
        nP++;
      }
      break;

    case GLMM::BERNOULLI_LOGIT:
      for (i = 0; i < *I; i++){
        for (j = 0; j < *nP; j++){
          *dYP    = 0;
          *meanYP = AK_Basic::invlogit_AK(*etaP);

          Y_dP++;
          dYP++;
          meanYP++;
          etaP++;
        }
        nP++;
      }
      break;

    case GLMM::POISSON_LOG:
      sum_dY_iP = sum_dY_i;
      for (i = 0; i < *I; i++){
        for (j = 0; j < *nP; j++){
          *dYP    = lgamma1p(double(*Y_dP));    /* = log(Gamma(1 + Y_d)) = log(Y_d!) */
          *meanYP = AK_Basic::exp_AK(*etaP);          
          *sum_dY_iP += *dYP;  

          Y_dP++;
          dYP++;
          meanYP++;
          etaP++;
        }
        sum_dY_iP++;
        nP++;
      }     
      break;

    default:
      *err = 1;
      //Rf_error("GLMM::dY_meanY: Unimplemented distributional type.\n", *distP);
      Rf_error("GLMM::dY_meanY: Unimplemented distributional type.\n");                 /* replaced the previous row on 08/12/2023 */
    }

    distP++;
  }  

  *sum_dY = AK_Basic::sum(sum_dY_i, *I);

  return;
}

}    // end of namespace GLMM
