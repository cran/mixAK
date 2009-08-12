//
//  PURPOSE:   Implementation of methods declared in GLMM_updateVars_eps.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/07/2009
//
// ======================================================================
//
#include "GLMM_updateVars_eps.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateVars_eps                                                                      *****/
/***** ***************************************************************************************** *****/
void
updateVars_eps(double* sigma,  
               const double* Y,     const double* eta,  
               const int* R,        const int* I,            const int* n,  const int* N_s,
               const double* zeta,  const double* gammaInv)
{
  static int s, i, j;
  static double resid, shape, scale;
  static double *sigmaP;
  static const double *zetaP, *gammaInvP;
  static const double *YP, *etaP;
  static const int *nP, *N_sP;

  sigmaP    = sigma;
  zetaP     = zeta;
  gammaInvP = gammaInv;
  YP        = Y;
  etaP      = eta;
  nP        = n;
  N_sP      = N_s;
  for (s = 0; s < *R; s++){

    /** Compute sum[i]sum[j](y[s,i,j] - eta[s,i,j])^2 **/
    scale = 0.0;
    for (i = 0; i < *I; i++){
      for (j = 0; j < *nP; j++){    /** *nP can be zero as well, does not matter here **/
        resid = *YP - *etaP;
        scale += resid * resid; 
        YP++;
        etaP++;
      }
      nP++;
    }

    /** Scale and shape of the full conditional gamma distribution **/
    scale = 1 / (0.5 * (*gammaInvP + scale));
    shape = 0.5 * (*zetaP + *N_sP);

    /** Sample new inverse variance and transform it to std. deviation **/
    *sigmaP = rgamma(shape, scale);
    *sigmaP = 1 / sqrt(*sigmaP);

    sigmaP++;
    zetaP++;
    gammaInvP++;
    N_sP++;
  }

  return;
}


}    /*** end of namespace GLMM ***/
