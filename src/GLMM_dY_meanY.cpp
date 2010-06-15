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
         double* meanY,
         int*    err,
         const double* Y_c,
         const int*    Y_d,
         const double* eta,
         const int*    dist,
         const int*    N_s,
         const int*    R_c,
         const int*    R_d)
{
  int s, i;

  const double *Y_cP  = Y_c;
  const int    *Y_dP  = Y_d;
  const double *etaP  = eta;
  const int    *distP = dist;
  const int    *N_sP  = N_s;

  double *dYP    = dY;
  double *meanYP = meanY;

  for (s = 0; s < (*R_c + *R_d); s++){
    switch (*distP){
    case GLMM::GAUSS_IDENTITY:
      for (i = 0; i < *N_sP; i++){
        *dYP    = 0;
        *meanYP = *etaP;

        Y_cP++;
        dYP++;
        meanYP++;
        etaP++;
      }
      break;

    case GLMM::BERNOULLI_LOGIT:
      for (i = 0; i < *N_sP; i++){
        *dYP    = 0;
        *meanYP = AK_Basic::invlogit_AK(*etaP);

        Y_dP++;
        dYP++;
        meanYP++;
        etaP++;
      }
      break;

    case GLMM::POISSON_LOG:
      for (i = 0; i < *N_sP; i++){
        *dYP    = lgamma1p(double(*Y_dP));    /* = log(Gamma(1 + Y_d)) = log(Y_d!) */
        *meanYP = AK_Basic::exp_AK(*etaP);

        Y_dP++;
        dYP++;
        meanYP++;
        etaP++;
      }     
      break;

    default:
      *err = 1;
      error("GLMM::dY_meanY: Unimplemented distributional type.\n", *distP);
    }

    N_sP++;
    distP++;
  }  


  return;
}

}    // end of namespace GLMM
