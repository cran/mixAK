//
//  PURPOSE:   Implementation of methods declared in Dist_TmixMVN.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/11/2008
//
// ====================================================================================================
//
#include "Dist_TmixMVN.h"

namespace Dist{

/***** *********************************************************************************** *****/
/***** Dist::rTmixMVN1                                                                     *****/
/***** *********************************************************************************** *****/
void
rTmixMVN1(double* x,  
          const int* K,        const double* cumw,  
          const double* beta,  const double* sigmaR2,  
          const double* a,     const double* b,        const int* trunc,  
          const int* p,        const int* p_p)
{
  static int r;

  /*** Sample the component ***/
  Dist::rDiscrete2(&r, cumw, K);

  /*** Sample the value from the component ***/
  Dist::rTMVN1(x, beta + r*(*p_p), sigmaR2 + r*(*p), a, b, trunc, p);

  return;
}

}  /** end of namespace Dist **/
