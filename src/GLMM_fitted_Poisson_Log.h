//
//  PURPOSE:   GLM, compute fitted values (Poisson expectations) for
//             GLM with Poisson distribution and log link 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   21/10/2009
//
//  FUNCTIONS:  
//     *   21/10/2009:  GLMM::fitted_Poisson_Log, PROTOTYPE 1
//     *   21/10/2009:  GLMM::fitted_Poisson_Log, PROTOTYPE 2
//
// ======================================================================
//
#ifndef _GLMM_FITTED_POISSON_LOG_H_
#define _GLMM_FITTED_POISSON_LOG_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Poisson_Log                                                                  *****/
/*****   PROTOTYPE 1                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  fitted[nobs]:   INPUT:  whatsever
//                 OUTPUT:  computed fitted values = exp(eta_fixed + eta_random)
//
//  eta_fixed[nobs]:   first part of the linear predictor
//
//  eta_random[nobs]:  second part of the linear predictor
//
//  nobs[1]:           number of observations
//
void
fitted_Poisson_Log(double* fitted,
                   const double* eta_fixed,  const double* eta_random, 
                   const int* nobs);


/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Poisson_Log                                                                  *****/
/*****   PROTOTYPE 2                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  fitted[sum(n)]:   INPUT:  whatsever
//                   OUTPUT:  computed fitted values = exp(eta_fixed + eta_random)
//
//  eta_fixed[sum(n)]:   first part of the linear predictor
//
//  eta_random[sum(n)]:  second part of the linear predictor
//
//  I[1]:                number of clusters
//
//  n[I]:                number of observations within each cluster
//
void
fitted_Poisson_Log(double* fitted,
                   const double* eta_fixed,  const double* eta_random, 
                   const int* I,             const int* n);

}

#endif
