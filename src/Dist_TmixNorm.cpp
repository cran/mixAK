//
//  PURPOSE:   Implementation of methods declared in Dist_TmixNorm.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/11/2008
//
// ====================================================================================================
//
#include "Dist_TmixNorm.h"

namespace Dist{

/***** ************************************************************************** *****/
/***** Dist::rTmixNorm1                                                           *****/
/***** ************************************************************************** *****/
void
rTmixNorm1(double* x,
           const int* K,     const double* cumw,  const double* mu,  const double* sigma,
           const double* a,  const double* b,     const int* trunc)
{
  static int r;

  /*** Sample the component ***/
  Dist::rDiscrete2(&r, cumw, K);

  /*** Sample the value from the component ***/
  Dist::rTNorm1(x, mu + r, sigma + r, a, b, trunc);  
  
  return;
}

}  /** end of namespace Dist **/
