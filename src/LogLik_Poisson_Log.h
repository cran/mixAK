//
//  PURPOSE:   Log-likelihood, score vector, information matrix (in this case, observed = expected)
//             for Poisson GLM with log link
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
//  PART OF CODE TAKEN FROM:   R package glmmAK, ll_poisson.{cpp,h}
//
//  FUNCTIONS:  
//     *   19/10/2009:  LogLik::Poisson_Log (PROTOTYPE 1)
//     *   19/10/2009:  LogLik::Poisson_Log (PROTOTYPE 2)
//
// =================================================================================
//
#ifndef _LOGLIK_POISSON_LOG_H_
#define _LOGLIK_POISSON_LOG_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace LogLik{

/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log (PROTOTYPE 1)                                                         *****/
/***** ***************************************************************************************** *****/
//
// This prototype
//   * updates the linear predictor
//   * updates the Poisson means
//   * computes log-likelihood, score and information matrix
//
//  In the case that we have random effects b = shift + scale * b^*
//  and theta = b
//  then the score vector and the information matrix are computed as derivatives
//  of the log-likelihood with respect to b^*
//
//  ll[1]:                    INPUT:  whatsever
//                           OUTPUT:  computed value of the log-likelihood
//
//  U[p + Intcpt]:            INPUT:  whatsever
//                           OUTPUT:  computed score vector
//
//  I[LT(p + Intcpt)]:        INPUT:  whatsever
//                           OUTPUT:  computed information matrix (lower triangle in COLUMN major order)
//                                    in this case, observed information matrix = expected information matrix
//
//  eta[n]:                   INPUT:  whatsever
//                           OUTPUT:  current value of the linear predictor related to the regression coefficients
//                                    stored in theta
//                                    eta[i] = t(x[,i]) %*% theta
//
//  lambda[n]:                INPUT:  whatsever
//                           OUTPUT:  current value of the Poisson mean
//                                    = exp(eta + offset)
//
//  offset[n]:                 offset vector
//                             In MCMC applications, this will typically be the part of the linear predictor
//                             not related to theta
//  
//  theta[p + Intcpt]:         regression coefficients
//                             mostly:  fixed effects  alpha
//                             or       random effects b = shift + scale * b^*
//                  
//  y[n]:                      observed counts
//
//  log_y_factor[n]:           log(y!) = lgamma(y + 1) for each observation
//
//  scale[p + Intcpt]:         scale factor used to get theta
//                             mostly:   for fixed effects:  scale = (1, 1, ..., 1)
//                                       for random effects: scale = scale vector used for random effects
//
//  x[p, n]:                   covariates (without intercept) for each observation
//
//  SxxS[LT(p + Intcpt), n]:   lower triangles of matrices t(S) %*% x %*% t(x) %*% S for each observation
//                             where S is a diagonal scale matrix, that is (in standard matrix notation)
//                             SxxS[j,k; n] = scale[j] * x[j, n] * x[k, n] * scale[k]
//
//  n[1]:                      number of observations
//
//  p[1]:                      number of covariates (intercept excluding)
//
//  Intcpt[1]:                 0/1 indicating whether intercept is included in the model
//
void
Poisson_Log1(double* ll,
             double* U,
             double* I,
             double* eta,
             double* lambda,
             const double* offset,
             const double* theta,
             const int* y,
             const double* log_y_factor,
             const double* scale,
             const double* x,
             const double* SxxS,
             const int* n,
             const int* p,
             const int* Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log (PROTOTYPE 2)                                                         *****/
/***** ***************************************************************************************** *****/
//
// This prototype only computes log-likelihood, score and information matrix
// from supplied eta, offset, lambda
//
// It is assumed that log(lambda) = eta + offset
//
void
Poisson_Log2(double* ll,
             double* U,
             double* I,
             const double* eta,
             const double* offset,
             const double* lambda,
             const int* y,
             const double* log_y_factor,
             const double* scale,
             const double* x,
             const double* SxxS,
             const int* n,
             const int* p,
             const int* Intcpt);

}

#endif
