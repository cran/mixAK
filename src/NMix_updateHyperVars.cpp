//
//  PURPOSE:   Implementation of methods declared in NMix_updateHyperVars.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2007
//
// ======================================================================
//
#include "NMix_updateHyperVars.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateHyperVars                                                                     *****/
/***** ***************************************************************************************** *****/
void
updateHyperVars(double* gammaInv,    double* XiInv,    double* log_sqrt_detXiInv,  double* dwork,
                const double* Q,     const int* K,     const int* p,
                const double* zeta,  const double* g,  const double* h)
{
  static int j, l;
  static double shape, scale;
  static double *gammaInvP, *XiInvP, *sumQllP;
  static const double *QP, *gP, *hP;

  /*** Fill dwork (sumQll) by zeros ***/
  sumQllP = dwork;
  for (l = 0; l < *p; l++){
    *sumQllP = 0.0;
    sumQllP++;
  }
  
  /*** Compute sum(Q_j[l,l]) (j=0, ..., K-1) ***/
  QP = Q;
  for (j = 0; j < *K; j++){
    sumQllP = dwork;
    for (l = 0; l < *p; l++){
      *sumQllP += *QP;
      sumQllP++;
      QP += *p - l; 
    }
  }

  /*** Sample new gamma^{-1}'s      ***/
  /*** Compute log_sqrt_detXiInv    ***/
  gammaInvP = gammaInv;
  XiInvP    = XiInv;
  sumQllP   = dwork; 
  gP        = g;
  hP        = h;
  *log_sqrt_detXiInv = 0.0;
  for (l = *p; l > 0; l--){
    shape = *gP + 0.5 * *K * *zeta;
    scale = 1/(*hP + 0.5 * *sumQllP);
    *gammaInvP = rgamma(shape, scale);

    *XiInvP = *gammaInvP;
    *log_sqrt_detXiInv += AK_Basic::log_AK(*gammaInvP);

    sumQllP++;
    gammaInvP++;
    XiInvP += l;
    gP++;
    hP++;
  }
  *log_sqrt_detXiInv *= 0.5;

  return;
}

}    /*** end of namespace NMix ***/

