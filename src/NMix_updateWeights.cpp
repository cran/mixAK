//
//  PURPOSE:   Implementation of methods declared in NMix_updateWeights.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007
//
// ======================================================================
//
#include "NMix_updateWeights.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateWeights                                                                       *****/
/***** ***************************************************************************************** *****/
void
updateWeights(double* w,           double* logw,  double* dwork,
              const int* mixN,     const int* K,  
              const double* delta, 
              const int* mixNxw,   const int* nxw)
{
  static int j, ixw;
  static double *logwP, *dworkP;
  static double *wP; 
  static const int *mixNP;

  static double *wxP;

  /*** No covariates on mixture weights ***/
  /*** ================================ ***/
  if (*nxw == 1){

    /*** Compute parameters of the Dirichlet distribution in the full conditional ***/
    dworkP = dwork;
    mixNP  = mixN;
    for (j = 0; j < *K; j++){
      *dworkP = *delta + *mixNP;
      dworkP++;
      mixNP++;
    }

    /*** Sample new weights ***/
    Dist::rDirichlet(w, dwork, K);

    /*** Recalculate log(w) ***/
    wP    = w;
    logwP = logw;
    for (j = 0; j < *K; j++){
      *logwP = AK_Basic::log_AK(*wP);
      wP++;
      logwP++;
    }
  }

  /*** Factor covariate on mixture weights (nxw > 1) ***/
  /*** ============================================= ***/
  else{

    wxP   = w;
    logwP = logw;
    mixNP = mixNxw;

    for (ixw = 0; ixw < *nxw; ixw++){

      /*** Compute parameters of the Dirichlet distribution in the full conditional ***/
      dworkP = dwork;
      for (j = 0; j < *K; j++){
        *dworkP = *delta + *mixNP;
        dworkP++;
        mixNP++;
      }

      /*** Sample new weights ***/
      Dist::rDirichlet(wxP, dwork, K);

      /*** Recalculate log(w) ***/
      wP    = wxP;
      for (j = 0; j < *K; j++){
        *logwP = AK_Basic::log_AK(*wP);
        wP++;
        logwP++;
      }

      /*** Move pointer to next weight set ***/
      wxP = wP;
    }
  }

  return;
}

}    /*** end of namespace NMix ***/
