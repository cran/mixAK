//
//  PURPOSE:   Implementation of methods declared in GLMM_updateHyperVars_eps.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/07/2009
//
// ======================================================================
//
#include "GLMM_updateHyperVars_eps.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateHyperVars_eps                                                                 *****/
/***** ***************************************************************************************** *****/
void
updateHyperVars_eps(double* gammaInv,  
                    const double* sigma,  const int* R,
                    const double* zeta,   const double* g,  const double* h)
{
  static int s;
  static double shape, scale;
  static double *gammaInvP;
  static const double *sigmaP, *zetaP, *gP, *hP;

  gammaInvP = gammaInv;
  sigmaP    = sigma;
  zetaP     = zeta;
  gP        = g;
  hP        = h;
  for (s = 0; s < *R; s++){
    shape = *gP + 0.5 * *zetaP;
    scale = 1 / (*hP + 0.5 * (1 / (*sigmaP * *sigmaP)));
    *gammaInvP = rgamma(shape, scale);

    gammaInvP++;
    sigmaP++;
    zetaP++;
    gP++;
    hP++;
  }

  return;
}

}    /*** end of namespace GLMM ***/
