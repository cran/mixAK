//
//  PURPOSE:   Implementation of methods declared in NMix_updateCensObs.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   27/11/2007
//
// ======================================================================
//
#include "NMix_updateCensObs.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateCensObs                                                                       *****/
/***** ***************************************************************************************** *****/
void
updateCensObs(double* y,         double* beta,      double* sigmaR2,      double* dwork,  int* err,          
              const double* y0,  const double* y1,  const int* censor,  
              const int* r,      const double* mu,  const double* Sigma,  
              const int* K,      const int* p,      const int* n)
{
  static int i, j, LTp, p_p;
  static double *yP, *sigmaP, *betaP, *sigmaR2P;
  static const double *y0P, *y1P, *SigmaP, *muP;
  static const int *censorP, *rP;

  if (*p == 1){

    /*** Compute mixture standard deviations ***/
    sigmaP = dwork;
    SigmaP = Sigma;
    for (j = 0; j < *K; j++){
      *sigmaP = sqrt(*SigmaP);
      sigmaP++;
      SigmaP++;
    }

    /*** Sample ***/
    yP      = y;
    y0P     = y0;
    y1P     = y1;
    censorP = censor;
    rP      = r;
    for (i = 0; i < *n; i++){
      muP    = mu + *rP;
      sigmaP = dwork + *rP;
      Dist::rTNorm1(yP, muP, sigmaP, y0P, y1P, censorP);
      yP++;
      y0P++;
      y1P++;
      censorP++;
      rP++;
    }
    return;
  }

  /***** ========== *p > 1 ========== *****/
  LTp = (*p * (*p + 1))/2;
  p_p = *p * *p;

  /*** Compute regressions to be able to get conditional distributions ***/
  betaP    = beta;
  sigmaR2P = sigmaR2;   
  muP      = mu;
  SigmaP   = Sigma;
  for (j = 0; j < *K; j++){
    Stat::BLA(betaP, sigmaR2P, dwork, err, muP, SigmaP, p);
    betaP    += p_p;
    sigmaR2P += *p;
    muP      += *p;
    SigmaP   += LTp;
  }

  /*** Sample ***/
  yP      = y;
  y0P     = y0;
  y1P     = y1;
  censorP = censor;
  rP      = r;
  for (i = 0; i < *n; i++){
    betaP    = beta + *rP * p_p;
    sigmaR2P = sigmaR2 + *rP * *p;
    Dist::rTMVN1(yP, betaP, sigmaR2P, y0P, y1P, censorP, p);
    yP      += *p;
    y0P     += *p;
    y1P     += *p;
    censorP += *p;
    rP++;
  }

  return;
}

}    /*** end of namespace NMix ***/


