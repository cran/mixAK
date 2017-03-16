//
//  PURPOSE:   Normal mixture model, the third part of the log-Jacobian
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/01/2008
//
//  FUNCTIONS:  
//     *  ??/01/2008:   NMix::logJac_part3
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_LOG_JAC_PART_THREE_H_
#define _NMIX_RJMCMC_LOG_JAC_PART_THREE_H_

#include <R.h>
#include <Rmath.h>        // added on 15/03/2017

#include "AK_Basic.h"

namespace NMix{

/***** ****************************************************************************************************************** *****/
/***** NMix::RJMCMC_logJac_part3:  The third part of the log-Jacobian                                                     *****/
/***** ****************************************************************************************************************** *****/
//
// logJac3[1]:      Computed third part of the log-Jacobian
//
// Lambdastar[p]:   Eigenvalues (in descending order) of the old splitted or new combined component
//
// Vstar[p, p]:     Eigenvectors (in columns) of the old splitted or new combined component
//
// p[1]:            Dimension
//
void
RJMCMC_logJac_part3(double* logJac3,  const double* Lambdastar,  const double* Vstar,  const double* P,  const int* p);

}    /*** end of namespace NMix ***/

#endif
