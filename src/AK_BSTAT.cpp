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


/***** *************************************************************************************************** *****/
/***** AK_BSTAT::shiftScale                                                                                *****/
/***** *************************************************************************************************** *****/
void
shiftScale(double* yscaled,  const double* y,  const double* shift,  const double* scale,  const int* n,  const int* p)
{
  static const double *yP;
  static const double *shiftP;
  static const double *scaleP;
  static double *yscaledP;
  static int i, j;

  yP       = y;
  yscaledP = yscaled;
  for (i = 0; i < *n; i++){
    shiftP = shift;
    scaleP = scale;
    for (j = 0; j < *p; j++){
      *yscaledP = (*yP - *shiftP) / *scaleP;
      yscaledP++;
      yP++;
      shiftP++;
      scaleP++;
    }
  }
  
  return;
}


/***** *************************************************************************************************** *****/
/***** AK_BSTAT::inv_shiftScale                                                                            *****/
/***** *************************************************************************************************** *****/
void
inv_shiftScale(double* y,  const double* yscaled,  const double* shift,  const double* scale,  const int* n,  const int* p)
{
  static const double *yscaledP;
  static const double *shiftP;
  static const double *scaleP;
  static double *yP;
  static int i, j;

  yP       = y;
  yscaledP = yscaled;
  for (i = 0; i < *n; i++){
    shiftP = shift;
    scaleP = scale;
    for (j = 0; j < *p; j++){
      *yP = *shiftP + *scaleP * *yscaledP;
      yscaledP++;
      yP++;
      shiftP++;
      scaleP++;
    }
  }
  
  return;
}


}  /*** end of namespace AK_BSTAT ***/
