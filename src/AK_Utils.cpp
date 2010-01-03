//
//  PURPOSE:   Implementation of methods declared in AK_Utils.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/11/2007
//
// ======================================================================
//
#include "AK_Utils.h"

namespace AK_Utils{

/***** ***************************************************************************************** *****/
/***** AK_Utils::R_rsort_desc                                                                    *****/
/***** ***************************************************************************************** *****/
void
R_rsort_desc(double* a,  const int& n)
{
  static int i;
  static double *aP;

  /** a = -a **/
  aP = a;
  for (i = 0; i < n; i++){
    *aP *= (-1);
    aP++;
  }

  /** Sort -a in ascending order **/
  R_rsort(a, n);

  /** a = -a **/
  aP = a;
  for (i = 0; i < n; i++){
    *aP *= (-1);
    aP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** AK_Utils::cum_Pr2Pr                                                                       *****/
/***** ***************************************************************************************** *****/
void
cum_Pr2Pr(double* Pr,  const double* cum_Pr,
          const int* K,  const int *I)
{
  static int i, j;
  static double *PrP;
  static const double *cum_PrP;
  static const double *sum_PrP;

  PrP     = Pr;
  cum_PrP = cum_Pr;
  for (i = 0; i < *I; i++){
    sum_PrP = cum_PrP + (*K - 1);

    /** j = 0 **/
    *PrP = *cum_PrP / *sum_PrP;
    PrP++;
    cum_PrP++;

    /** j = 1, ..., K-1 **/
    for (j = 1; j < *K; j++){
      *PrP = (*cum_PrP - *(cum_PrP - 1)) / *sum_PrP;
      PrP++;
      cum_PrP++;
    }
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** AK_Utils::printIterInfo                                                                   *****/
/***** ***************************************************************************************** *****/
void
printIterInfo(int& writeAll,  int& backs,  const int& iter,  const int& nwrite,  const int& lastIter)
{
  static int i;

  if (!(iter % nwrite) || iter == lastIter){
    writeAll = 1;
    for (i = 0; i < backs; i++) Rprintf((char*)("\b"));
    Rprintf((char*)("%d"), iter);
    backs = int(log10(double(iter))) + 1;
  }

  return;
}

}    /*** end of namespace AK_Utils ***/
