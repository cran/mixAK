//
//  PURPOSE:   Basic statistics
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007 as AK_Utils.h
//             14/11/2007
//
//  FUNCTIONS:
//     * yBar_s  05/11/2007:   Compute sample means and standard deviations from a sample
//                             from a (multivariate) distribution. For both, ML etimates are provided,
//                             that is, we always divide by "n".
//
// ======================================================================
//
#ifndef _AK_BASIC_STATISTICS_H_
#define _AK_BASIC_STATISTICS_H_

#include <R.h>

namespace AK_BSTAT{

/***** *************************************************************************************************** *****/
/***** AK_BSTAT::yBar_s                                                                                    *****/
/***** *************************************************************************************************** *****/
// yBar[dimy[1]]         Computed sample means.
//
// ySD[dimy[1]]          Computed sample standard deviations.
//
// y[dimy[0], dimy[1]]   Matrix (in column major order) of data.
//
// dimy[2]               Dimension of the data
//                       * dimy[0] = number of observations
//                       * dimy[1] = dimension of the distribution
//
void
yBar_s(double* yBar,  double* ySD,  const double* y,  const int* dimy);


}  /*** end of namespace AK_BSTAT ***/

#endif
