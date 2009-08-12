//
//  PURPOSE:   (Multivariate) GLMM, compute values of linear predictors from scratch
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/07/2009
//
//  FUNCTIONS:  
//     *   12/07/2009:  GLMM::linear_predictors  
//     *   07/08/2009:  GLMM::linear_predictor_fixed
//     *   07/08/2009:  GLMM::linear_predictor_random
//     *   07/08/2009:  GLMM::linear_predictor_zs
//
// ======================================================================
//
#ifndef _GLMM_LINEAR_PREDICTORS_H_
#define _GLMM_LINEAR_PREDICTORS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors                                                                   *****/
/***** ***************************************************************************************** *****/
//
//  eta_fixed[sum(n)]:      INPUT:  whatsever
//                         OUTPUT:  values of x[s,i,j]'beta (+ intercept)
//                                  filled by zeros if there are no betas/intercepts                               
//  
//  eta_random[sum(n)]:     INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'b[s,i] (+ random intercept)
//                                  filled by zeros if there are no random effects
//
//  eta[sum(n)]:            INPUT:  whatsever
//                         OUTPUT:  values of eta_fixed + eta_random
//
//  eta_zs[sum(n)]:         INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'shift_b[s]
//
//  N_s[R]:                 INPUT:  whatsever
//                         OUTPUT:  N_s[s] = number of observations for response s
//                                  sum(N_s) = sum(n)
//
//  OTHER VARIABLES:  see GLMM_MCMC.{h,cpp}
//
void
linear_predictors(double* eta_fixed,  double* eta_random,      double* eta,      double* eta_zs,         int* N_s,
                  const double* X,    const double* beta,      const double* Z,  const double* b,        const double* shift_b,
                  const int* p,       const int* fixedIntcpt,  const int* q,     const int* randIntcpt,  
                  const int* n,       const int* R,            const int* I,     const int* dim_b,       const int* cumq_ri);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_fixed                                                              *****/
/***** ***************************************************************************************** *****/
//
//  eta_fixed[sum(n), I]:   INPUT:  whatsever
//                         OUTPUT:  values of x[s,i,j]'beta (+ intercept)
//                                  filled by zeros if there are no betas/intercepts
//  
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictors2 function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_fixed(double* eta_fixed,
                       const double* X,    const double* beta,
                       const int* p,       const int* fixedIntcpt,
                       const int* n,       const int* R,            const int* I);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_random                                                             *****/
/***** ***************************************************************************************** *****/
//
//  eta_random[sum(n), I]:  INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'b[s,i] (+ random intercept)
//                                  filled by zeros if there are no random effects
//
// 
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictors2 function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_random(double* eta_random,
                        const double* Z,     const double* b,
                        const int* q,        const int* randIntcpt,  
                        const int* n,        const int* R,            const int* I,     const int* dim_b,       const int* cumq_ri);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_zs                                                                 *****/
/***** ***************************************************************************************** *****/
//  eta_zs[sum(n), I]:      INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'shift_b[s]
//
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictors2 function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_zs(double* eta_zs, 
                    const double* Z,  const double* shift_b,
                    const int* q,     const int* randIntcpt,  
                    const int* n,     const int* R,           const int* I,        const int* dim_b,       const int* cumq_ri);

}

#endif
