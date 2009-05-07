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
