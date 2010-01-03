//
//  PURPOSE:   Log-likelihood, score vector, information matrix (in this case, observed = expected)
//             for Gaussian GLM with identity link (i.e., standard normal linear model)
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/10/2009
//
//  FUNCTIONS:  
//     *   19/10/2009:  LogLik::Gauss_Identity (PROTOTYPE 1)
//     *   03/11/2009:  LogLik::Gauss_Identity (PROTOTYPE 3)
//     *   27/10/2009:  LogLik::Gauss_Identity (PROTOTYPE 4)
//
// =================================================================================
//
#ifndef _LOGLIK_GAUSS_IDENTITY_H_
#define _LOGLIK_GAUSS_IDENTITY_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace LogLik{

/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity (PROTOTYPE 1)                                                      *****/
/***** ***************************************************************************************** *****/
//
// This prototype
//   * updates the linear predictor
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
//  mu[n]:                    INPUT:  whatsever
//                           OUTPUT:  current value of the normal mean
//                                    = eta + offset
//
//  offset[n]:                 offset vector
//                             In MCMC applications, this will typically be the part of the linear predictor
//                             not related to theta
//  
//  theta[p + Intcpt]:         regression coefficients
//                             mostly:  fixed effects  alpha
//                             or       random effects b = shift + scale * b^*
//                  
//  y[n]:                      observed outcome
//
//  sigma[1]:                  residual standard deviation
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
Gauss_Identity1(double* ll,
                double* U,
                double* I,
                double* eta,
                double* mu,
                const double* offset,
                const double* theta,
                const double* y,
                const double* sigma,
                const double* scale,
                const double* x,
                const double* SxxS,
                const int* n,
                const int* p,
                const int* Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity (PROTOTYPE 3)                                                      *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes the log-likelihood and updates the linear predictor eta
//
void
Gauss_Identity3(double* ll,
                double* eta,
                const double* offset,
                const double* theta,
                const double* y,
                const double* sigma,
                const double* x,
                const int* n,
                const int* p,
                const int* Intcpt);


/***** ***************************************************************************************** *****/
/***** LogLik::Gauss_Identity (PROTOTYPE 4)                                                      *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes only the log-likelihood from supplied y and eta
//
void
Gauss_Identity4(double* ll,
                const double* eta,
                const double* offset,
                const double* y,
                const double* sigma,
                const int* n);

}

#endif
