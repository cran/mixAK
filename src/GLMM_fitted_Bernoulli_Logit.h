//
//  PURPOSE:   GLM, compute fitted values (probabilities of success) for
//             GLM with Bernoulli (alternative) distribution and logit link 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   21/10/2009
//
//  FUNCTIONS:  
//     *   21/10/2009:  GLMM::fitted_Bernoulli_Logit, PROTOTYPE 1
//     *   21/10/2009:  GLMM::fitted_Bernoulli_Logit, PROTOTYPE 2
//
// ======================================================================
//
#ifndef _GLMM_FITTED_BERNOULLI_LOGIT_H_
#define _GLMM_FITTED_BERNOULLI_LOGIT_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Bernoulli_Logit                                                              *****/
/*****   PROTOTYPE 1                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  fitted[nobs]:   INPUT:  whatsever
//                 OUTPUT:  computed fitted values = invlogit(eta_fixed + eta_random)
//
//  eta_fixed[nobs]:   first part of the linear predictor
//
//  eta_random[nobs]:  second part of the linear predictor
//
//  nobs[1]:           number of observations
//
void
fitted_Bernoulli_Logit(double* fitted,
                       const double* eta_fixed,  const double* eta_random, 
                       const int* nobs);


/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Bernoulli_Logit                                                              *****/
/*****   PROTOTYPE 2                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  fitted[sum(n)]:   INPUT:  whatsever
//                   OUTPUT:  computed fitted values = invlogit(eta_fixed + eta_random)
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
fitted_Bernoulli_Logit(double* fitted,
                       const double* eta_fixed,  const double* eta_random, 
                       const int* I,             const int* n);

}

#endif
