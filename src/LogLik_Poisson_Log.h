//
//  PURPOSE:   Log-likelihood, score vector, information matrix (in this case, observed = expected)
//             for Poisson GLM with log link
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
//  LOG:       06/12/2011:  R_NegInf replaced by log(0) = AK_Basic::_LOG_ZERO0 when calculating log-likelihood
//
//  PART OF CODE TAKEN FROM:   R package glmmAK, ll_poisson.{cpp,h}
//
//  FUNCTIONS:  
//     *   14/04/2010:  Poisson_Log1
//     *   09/04/2010:  Poisson_Log_sqrt_w_phi_stres1
//     *   09/04/2010:  Poisson_Log_sqrt_w_phi_stres2
//     *   19/10/2009:  LogLik::Poisson_LogUI1
//     *   19/10/2009:  LogLik::Poisson_LogUI1
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
//
//  ll[1]:                 INPUT:  whatsever
//                        OUTPUT:  computed value of the log-likelihood
//
//  sqrt_w_phi[n]:         INPUT:  whatsever
//                        OUTPUT:  1 * sqrt(lambda)
//             
//  stres[n]:              INPUT:  whatsever
//                        OUTPUT:  (y[i] - lambda[i]) / sqrt(lambda[i])
//
//  eta[n]:                
//
//  lambda[n]:
//
//  offset[n]:
//
//  theta[p + Intcpt]:
//
//  sqrt_phi[0]:           not used here (it is constantly equal to 1)
//
//  y[n]:
//
//  log_y_factor[n]:       log(y!) = lgamma(y + 1) for each observation                      
//
//  x[p, n]:
//
//  n[1]:
//
//  p[1]:
//
//  Intcpt[1]:
//
/***** ***************************************************************************************** *****/

/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log1                                                                      *****/
/***** ***************************************************************************************** *****/
//
// This version computes only the log-likelihood from offset, x, theta.
//
void
Poisson_Log1(double* ll,
             const double* offset,
             const double* theta,
             const double* sqrt_phi,
             const int*    y,
             const double* log_y_factor,
             const double* x,
             const int*    n,
             const int*    p,
             const int*    Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi1                                                           *****/
/***** ***************************************************************************************** *****/
//
//  This version computes log-likelihood and
//  sqrt_w_phi = phi^{-1} * sqrt(var(y | eta)) = 1 * sqrt(lambda).
//  
void
Poisson_Log_sqrt_w_phi1(double* ll,
                        double* sqrt_w_phi,
                        const double* offset,
                        const double* theta,
                        const double* sqrt_phi,
                        const int*    y,
                        const double* log_y_factor,
                        const double* x,
                        const int*    n,
                        const int*    p,
                        const int*    Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi_stres1                                                     *****/
/***** ***************************************************************************************** *****/
//
// This version computes: 1) eta from x and theta
//                        2) pi from eta and offset
//                        3) ll (log-likelihood)
//                        4) sqrt_w_phi = phi^{-1} * sqrt(var(y | eta)) = 1 * sqrt(lambda)
//                        5) stres = (y - lambda) / sqrt(var(y | eta))
//
void
Poisson_Log_sqrt_w_phi_stres1(double* ll,
                              double* sqrt_w_phi,
                              double* stres,
                              double* eta,
                              double* lambda,                       
                              const double* offset,
                              const double* theta,
                              const double* sqrt_phi,
                              const int*    y,
                              const double* log_y_factor,
                              const double* x,
                              const int*    n,
                              const int*    p,
                              const int*    Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_Log_sqrt_w_phi_stres2                                                     *****/
/***** ***************************************************************************************** *****/
//
// This version computes: 1) ll (log-likelihood)
//                        2) sqrt_w_phi = phi^{-1} * sqrt(var(y | eta)) = 1 * sqrt(lambda)
//                        3) stres = (y - lambda) / sqrt(var(y | eta))
//
// For this version: eta, offset, sqrt_phi can be NULL
//         
void
Poisson_Log_sqrt_w_phi_stres2(double* ll,
                              double* sqrt_w_phi,
                              double* stres,
                              const double* eta,
                              const double* offset,
                              const double* lambda,
                              const double* sqrt_phi,
                              const int*    y,
                              const double* log_y_factor,
                              const int*    n);


/***** ***************************************************************************************** *****/
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
/***** ***************************************************************************************** *****/

/***** ***************************************************************************************** *****/
/***** LogLik::Poisson_LogUI1                                                                    *****/
/***** ***************************************************************************************** *****/
//
// This version
//   * updates the linear predictor
//   * updates the Poisson means
//   * computes log-likelihood, score and information matrix
//
void
Poisson_LogUI1(double* ll,
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
/***** LogLik::Poisson_LogUI2                                                                    *****/
/***** ***************************************************************************************** *****/
//
// This version only computes log-likelihood, score and information matrix
// from supplied eta, offset, lambda
//
// It is assumed that log(lambda) = eta + offset
//
void
Poisson_LogUI2(double* ll,
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
