//
//  PURPOSE:   (Multivariate) GLMM.
//             Fill dY and meanY vectors (see GLMM_MCMC.cpp)
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   26/02/2010  (as file GLMM_dY_d_mean_Y_d.h by taking the code originally included in GLMM_MCMC.cpp)
//             12/04/2010  (as this file)
//
//  FUNCTIONS:  
//     * dY_meanY     12/04/2010
//
// ======================================================================
//
#ifndef _GLMM_DY_MEANY_H_
#define _GLMM_DY_MEANY_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "GLMM.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM:dY_meanY                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  dY[]               INPUT:  whatsever
//                    OUTPUT:  see GLMM_MCMC.cpp
//
//  meanY[]            INPUT:  whatsever
//                    OUTPUT:  see GLMM_MCMC.cpp
// 
//  err[1]            error flag
//
//  Y_c[]             continuous response variables
//
//  Y_d[]             discrete response variables
//
//  eta[]             values of linear predictors for both continuous and discrete responses
//
//  N_s[R_c + R_d]    numbers of observations for each response
//
//  R_c[1]            number of continuous responses
//
//  R_d[1]            number of discrete responses
//
/***** ***************************************************************************************** *****/
void
dY_meanY(double* dY,
         double* meanY,
         int*    err,
         const double* Y_c,
         const int*    Y_d,
         const double* eta,
         const int*    dist,
         const int*    N_s,
         const int*    R_c,
         const int*    R_d);

}    // end of namespace GLMM

#endif
