//
//  PURPOSE:   Implementation of methods declared in AK_BSTAT.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.cpp
//             14/11/2007
//
// ======================================================================
//
#include "AK_BSTAT.h"

namespace AK_BSTAT{

/***** ***************************************************************************************** *****/
/***** AK_BSTAT::yBar_s                                                                          *****/
/***** ***************************************************************************************** *****/
void
yBar_s(double* yBar,  double* ySD,  const double* y,  const int* dimy)
{
  static int i, j;
  static double tmp;
  static double *yBarP, *ySDP;
  static const double *yP;

  /*** Compute means ***/
  yBarP = yBar;
  yP = y;
  for (j = 0; j < dimy[1]; j++){
    *yBarP = 0.0;
    for (i=0; i < dimy[0]; i++){
      *yBarP += *yP;
      yP++;
    }
    *yBarP /= dimy[0];
    yBarP++;
  }

  /*** Compute standard deviations ***/
  yBarP = yBar;
  ySDP = ySD;
  yP = y;
  for (j=0; j<dimy[1]; j++){
    *ySDP = 0.0;
    for (i=0; i < dimy[0]; i++){
      tmp = *yP - *yBarP;
      *ySDP += tmp*tmp;
      yP++;
    }
    *ySDP /= dimy[0];
    *ySDP = sqrt(*ySDP);
    yBarP++;
    ySDP++;
  }

  return;
}

}  /*** end of namespace AK_BSTAT ***/
