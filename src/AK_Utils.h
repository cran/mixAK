//
//  PURPOSE:   Smaller utilities
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:    26/11/2007
//
//  FUNCTIONS:  
//       * R_rsort_desc   28/01/2008
//       * printIterInfo  26/11/2007
//
// ========================================================
//
#ifndef _AK_UTILS_H_
#define _AK_UTILS_H_

#include <R.h>

namespace AK_Utils{

/***** ***************************************************************************************** *****/
/***** AK_Utils::R_rsort_desc                                                                    *****/
/***** ***************************************************************************************** *****/
//
// Sort elements of a in DESCENDING order
// * main sorting is performed by R_rsort (R_ext/Utils.h included by R.h)
//
// a[n]   INPUT:   array to sort
//       OUTPUT:   sorted array
//
// n     length of a
//
void
R_rsort_desc(double* a,  const int& n);


/***** ***************************************************************************************** *****/
/***** AK_Utils::printIterInfo                                                                   *****/
/***** ***************************************************************************************** *****/
//
// Print information over performed MCMC on a screen
//
void
printIterInfo(int& writeAll,  int& backs,  const int& iter,  const int& nwrite,  const int& lastIter);

}  /*** end of namespace AK_Utils ***/

#endif


