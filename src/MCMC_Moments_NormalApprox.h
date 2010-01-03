//
//  PURPOSE:   Compute a canonical mean and a Cholesky decomposition of the precission matrix 
//             of the normal approximation based on one Newton-Raphson/Fisher scoring step 
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/10/2009
//
//  FUNCTIONS:  
//     *   26/10/2009:  MCMC:Moments_NormalApprox (PROTOTYPE 1)
//     *   27/10/2009:  MCMC:Moments_NormalApprox (PROTOTYPE 2)
//
// =================================================================================
//
#ifndef _MCMC_MOMENTS_NORMALAPPROX_H_
#define _MCMC_MOMENTS_NORMALAPPROX_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

namespace MCMC{

/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox (PROTOTYPE 1)                                                  *****/
/***** ***************************************************************************************** *****/
//
//  This prototype assumes that the prior covariance is diagonal and computes
//  * canonical mean of the normal approximation
//  * Cholesky decomposition of the precision matrix of the normal approximation
//  * log(|Q|^{1/2}), where Q is the precision matrix of the normal approximation
//
//  cmean[dim]:         INPUT:  score vector U
//                     OUTPUT:  computed canonical mean
//
//  Q[LT(dim)]:         INPUT:  information matrix I
//                     OUTPUT:  computed Cholesky decomposition of the precision matrix
//
//  log_sqrtdet_Q[1]:   INPUT:  whatsever
//                     OUTPUT:  log(|Q|^{1/2})
//
//  dwork[dim]:    working array
//
//  err[1]:        error flag
//                OUTPUT:  the same as input if no problems
//                         <> 0 if problems
//
//  theta[dim]:    current value of the parameter for which we construct 
//                 the normal approximation
//
//  Pprior[dim]:   vector of prior precisions
//                 It is assumed that the prior covariance matrix = diag(Pprior^{-1})
//
//  P_Mprior[dim]: vector with prior precision %*% prior mean
//
//  dim[1]:        dimension
//
//  caller:        name of the routine which has called this function 
//                 (used in error messages)
//
// *********************************************************************************************
void
Moments_NormalApprox(double* cmean, 
                     double* Q,
                     double* log_sqrtdet_Q,
                     double* dwork,
                     int* err,
                     const double* theta,
                     const double* Pprior,
                     const double* P_Mprior,
                     const int* dim,
                     const char* caller);


/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox (PROTOTYPE 2)                                                  *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes only canonical mean of the normal approximation
//
//  cmean[dim]:         INPUT:  score vector U
//                     OUTPUT:  computed canonical mean
//
//  dwork[dim]:    working array
//
//  theta[dim]:    current value of the parameter for which we construct 
//                 the normal approximation
//
//  Imat[LT(dim)]: information matrix
//
//  P_Mprior[dim]: vector with prior precision %*% prior mean
//
//  dim[1]:        dimension
//
// *********************************************************************************************
void
Moments_NormalApprox(double* cmean, 
                     double* dwork,                       
                     const double* theta,
                     const double* Imat,
                     const double* P_Mprior,
                     const int* dim);

}

#endif
