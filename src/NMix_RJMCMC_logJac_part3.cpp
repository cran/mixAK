//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMC_logJac_part3.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMC_logJac_part3.h"

namespace NMix {

/***** ****************************************************************************************************************** *****/
/***** NMix::RJMCMC_logJac_part3:  The third part of the log-Jacobian                                                     *****/
/***** ****************************************************************************************************************** *****/
void
RJMCMC_logJac_part3(double* logJac3,  const double* Lambdastar,  const double* Vstar,  const double* P,  const int* p)
{
  static double v22, v21;
  static const double *VP;

  switch (*p){
  case 1:
    *logJac3 = 0.0;
    return;

  case 2:             /*** Expression taken from Matlab code rj2logjacobian.m of I. Papageorgiou ***/
    VP = Vstar + 1;
    v21 = *VP;
    v22 = *(VP + 2);
    *logJac3 = M_LN2 + 2 * AK_Basic::log_AK(fabs(*P)) + AK_Basic::log_AK(v21*v21 + v22*v22) + AK_Basic::log_AK(fabs(v21));
    return;

  default:
    Rf_error("NMix::RJMCMC_logJac_part3 not (yet) implemented for dimension %d.\n", *p);
  }  
}


}    /*** end of namespace NMix ***/

