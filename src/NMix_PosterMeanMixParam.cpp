//
//  PURPOSE:   Implementation of methods declared in NMix_PosterMeanMixParam.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/02/2010
//
// ======================================================================
//
#include "NMix_PosterMeanMixParam.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::PosterMeanMixParam                                                                  *****/
/***** ***************************************************************************************** *****/
void
PosterMeanMixParam(double* pm_w,
                   double* pm_mu,
                   double* pm_Q,
                   double* pm_Sigma,
                   double* pm_Li,
                   const int*    K,
                   const double* chw,
                   const double* chmu,
                   const double* chQ,
                   const double* chSigma,
                   const double* chLi,
                   const int*    chorder,
                   const int*    p,
                   const int*    Mkeep)
{
  int i, j, l, l2;
  double *pm_wP, *pm_muP, *pm_QP, *pm_SigmaP, *pm_LiP;

  const int LTp = (*p * (*p + 1))/2;

  /*** Reset posterior means ***/
  AK_Basic::fillArray(pm_w,     0.0, *K);
  AK_Basic::fillArray(pm_mu,    0.0, *K * *p);
  AK_Basic::fillArray(pm_Q,     0.0, *K * LTp);
  AK_Basic::fillArray(pm_Sigma, 0.0, *K * LTp);
  AK_Basic::fillArray(pm_Li,    0.0, *K * LTp);

  /*** Initialize pointers ***/
  const double *chwP      = chw;
  const double *chmuP     = chmu;
  const double *chQP      = chQ;
  const double *chSigmaP  = chSigma;
  const double *chLiP     = chLi;
  const int    *chorderP  = chorder;  

  const double *chmuP2, *chQP2, *chSigmaP2, *chLiP2;  

  /*** Sums over sampled values***/
  for (i = 0; i < *Mkeep; i++){
    pm_wP     = pm_w;
    pm_muP    = pm_mu;
    pm_QP     = pm_Q;
    pm_SigmaP = pm_Sigma;
    pm_LiP    = pm_Li;

    for (j = 0; j < *K; j++){
      *pm_wP += chwP[*chorderP];
      pm_wP++;

      chmuP2    = chmuP    + (*chorderP * *p);
      chQP2     = chQP     + (*chorderP * LTp);
      chSigmaP2 = chSigmaP + (*chorderP * LTp);
      chLiP2    = chLiP    + (*chorderP * LTp);
      chorderP++;

      for (l2 = 0; l2 < *p; l2++){
        *pm_muP += *chmuP2;
        pm_muP++;
        chmuP2++;
        for (l = l2; l < *p; l++){
          *pm_QP     += *chQP2;
          *pm_SigmaP += *chSigmaP2;
          *pm_LiP    += *chLiP2;
          pm_QP++;
          pm_SigmaP++;
          pm_LiP++;
          chQP2++;
          chSigmaP2++;
          chLiP2++;
        }
      }
    }

    chwP     += *K;
    chmuP    += *p * *K;
    chQP     += LTp * *K;
    chSigmaP += LTp * *K;
    chLiP    += LTp * *K;
  }

  /*** Averages over sampled values ***/
  pm_wP     = pm_w;
  pm_muP    = pm_mu;
  pm_QP     = pm_Q;
  pm_SigmaP = pm_Sigma;
  pm_LiP    = pm_Li;
  for (j = 0; j < *K; j++){
    *pm_wP /= *Mkeep;
    pm_wP++;
    for (l2 = 0; l2 < *p; l2++){
      *pm_muP /= *Mkeep;
      pm_muP++;
      for (l = l2; l < *p; l++){
        *pm_QP     /= *Mkeep;
        *pm_SigmaP /= *Mkeep;
        *pm_LiP    /= *Mkeep;
        pm_QP++;
        pm_SigmaP++;
        pm_LiP++;
      }
    }
  }

  return;
}

}  // end of namespace NMix
