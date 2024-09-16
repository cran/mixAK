//
//  PURPOSE:   Multivariate normal distribution
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007
//             * some functions partially taken from GMRF.{h,cpp} of the glmmAK package 
//
//  FUNCTIONS:
//     * dMVN1  05/11/2007:   (Log-)density of the multivariate normal distribution N(mu, Q^{-1}),
//                             where the Cholesky decomposition Li of the precision matrix is given, i.e., Q = Li %*% t(Li).
//                             * partially taken from dGMRF function in GMRF.{h,cpp}[glmmAK]
//
//     * rMVN1  05/11/2007:    Random number generation from multivariate normal distribution N(mu, Q^{-1}),
//                             where the Cholesky decomposition Li of the precision matrix is given, i.e., Q = Li %*% t(Li).
//                             * function also returns log-density evaluated at sampled value
//                             * partially taken from rGMRF function in GMRF.{h,cpp}[glmmAK]
//
//     * rMVN2  12/11/2007:    Random number generation from multivariate normal distribution N(Q^{-1}*b, Q^{-1}),
//                             where the Cholesky decomposition Li of the precision matrix is given, i.e., Q = Li %*% t(Li),
//                             and a canonical mean b is given, i.e., the mean mu = Q^{-1}*b
//                             * function also returns log-density evaluated at sampled value
//
//     * rMVN3  10/11/2009:    Random number generation from multivariate normal distribution N(Q^{-1}*b, scale^{-1}*Q^{-1}),
//                             where the Cholesky decomposition Li of the precision matrix is given, i.e., Q = Li %*% t(Li),
//                             and a canonical mean b is given, i.e., the mean mu = Q^{-1}*b
//                             * function also returns log-density evaluated at sampled value
//
//     * rMVN4  12/04/2010:    Random number generation from multivariate normal distribution N(mu, scale^{-1}*Q^{-1}),
//                             where a factorization Li of the precision matrix is given, i.e., Q = Li %*% t(Li),
//                             and a mean mu is given
//                             * function also returns log-density evaluated at sampled value
//
//     * ldMVN1   03/12/2007:    Log-density of the multivariate normal distribution N(mu, Q^{-1})
//
//     * ldMVN2   03/12/2007:    Log-density of the multivariate normal distribution N(mu, Sigma)
//
//     * ldMVN3   10/11/2009:    Log-density of the multivariate normal distribution N(mu, scale^{-1}*Q^{-1})
//
//     * dMVN1_R  05/11/2007:    R wrapper for dMVN1 which allows to evaluate a density for several x points in a loop
//
//     * rMVN1_R  05/11/2007:    R wrapper for rMVN1 which allows to generate several x points in a loop
//                               * GetRNGstate() and PutRNGstate() is used inside!!!
//
//     * rMVN2_R  12/11/2007:    R wrapper for rMVN2 which allows to generate several x points in a loop
//                               * GetRNGstate() and PutRNGstate() is used inside!!!
//                     
// ==============================================================================================================================
//
#ifndef _DIST_MVN_H_
#define _DIST_MVN_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace Dist{


/***** ************************************************************************************************* *****/
/***** Dist::dMVN1                                                                                       *****/
/***** Dist::rMVN1                                                                                       *****/
/***** ************************************************************************************************* *****/
//
// (Log-)density of the multivariate normal distribution N(mu, Q^{-1})
//
// Random number generation from multivariate normal distribution N(mu, Q^{-1})
//
// log_dens[1]  Value of the (log-)density evaluated at x
// 
// work[nx]     Working array
//
// x[nx]        Point where the density is to be evaluated or generated random number
//
// unlog[1]     If not zero, then the value of the density is computed and not log-density
//
// mu[nx]       Mean of the normal distribution
//              * not needed when mu_nonZERO = 0 in which case it can be set to NULL
//
// Li[LT(nx)]   Cholesky decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// nx[1]          Dimension of the normal distribution
//
// mu_nonZERO[1]  If <> 0 then the argument 'mu' gives the mean of the normal distribution,
//                otherwise it is assumed that the mean of the normal distribution is zero
//
void
dMVN1(double* log_dens,  double* work,
      const double* x,   const int* unlog,
      const double* mu,  const double* Li,       const double* log_dets,
      const int* nx,     const int* mu_nonZERO);

void
rMVN1(double* x,         double* log_dens,
      const double *mu,  const double *Li,      const double *log_dets,           
      const int* nx,     const int* mu_nonZERO);


/***** ************************************************************************************************* *****/
/***** Dist::rMVN2                                                                                       *****/
/***** ************************************************************************************************* *****/
//
// Random number generation from multivariate normal distribution N(Q^{-1}*b, Q^{-1}),
//
// x[nx]        OUTPUT:  Generated random number
//
// mu[nx]       INPUT:   Canonical mean b of the normal distribution 
//              OUTPUT:  Computed mean of the normal distribution, i.e., mu = Q^{-1}*b
//
// log_dens[1]  Value of the (log-)density evaluated at x
// 
// work[nx]     Working array
//
// Li[LT(nx)]   Cholesky decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// nx[1]          Dimension of the normal distribution
//
void
rMVN2(double* x,         double* mu,              double* log_dens,
      const double *Li,  const double *log_dets,  const int* nx);


/***** ***************************************************************************************** *****/
/***** Dist::rMVN3                                                                               *****/
/***** ***************************************************************************************** *****/
//
// Random number generation from multivariate normal distribution N(Q^{-1}*b, scale^{-1}*Q^{-1}),
//
// x[nx]        OUTPUT:  Generated random number
//
// mu[nx]       INPUT:   Canonical mean b of the normal distribution 
//              OUTPUT:  Computed mean of the normal distribution, i.e., mu = Q^{-1}*b
//
// log_dens[1]  Value of the (log-)density evaluated at x
// 
// work[nx]     Working array
//
// Li[LT(nx)]   Cholesky decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2}) = sum(log(Li[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// sqrt_scale[1]          square root of the factor to scale the covariance matrix
//
// log_sqrt_scale[1]      log(sqrt(scale))
//
// nx[1]        Dimension of the normal distribution
//
void
rMVN3(double* x,         double* mu,              double* log_dens,
      const double *Li,  const double *log_dets,  const double *sqrt_scale,   const double *log_sqrt_scale,
      const int* nx);


/***** ***************************************************************************************** *****/
/***** Dist::rMVN4                                                                               *****/
/***** ***************************************************************************************** *****/
//
// Random number generation from multivariate normal distribution N(mu, scale^{-1}*Q^{-1}),
//
// x[nx]        OUTPUT:  Generated random number
//
// log_dens[1]  Value of the (log-)density evaluated at x
// 
// work[nx]     Working array
//
// mu[nx]       Mean of the normal distribution
//
// Li[LT(nx)]   Cholesky (or other) decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2}) = sum(log(Li[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// sqrt_scale[1]          square root of the factor to scale the covariance matrix
//
// log_sqrt_scale[1]      log(sqrt(scale))
//
// nx[1]        Dimension of the normal distribution
//
void
rMVN4(double* x,         double* log_dens,
      const double* mu,  const double *Li,  const double *log_dets,  const double *sqrt_scale,   const double *log_sqrt_scale,
      const int* nx);


/***** ***************************************************************************************** *****/
/***** Dist::ldMVN1                                                                              *****/
/***** Dist::ldMVN2                                                                              *****/
/***** ***************************************************************************************** *****/
//
// Log-density of the multivariate normal distribution N(mu, Q^{-1}) = N(mu, Sigma)
//
// log_dens[1]  Value of the log-density evaluated at x
// 
// work[nx]     Working array
//
// x[nx]        Point where the density is to be evaluated or generated random number
//
// mu[nx]       Mean of the normal distribution
//              * not needed when mu_nonZERO = 0 in which case it can be set to NULL
//
// Li[LT(nx)]   Cholesky decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// L[LT(nx)]    Cholesky decomposition of the covariance matrix Sigma (lower triangle only), i.e., Sigma = L %*% t(L)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j]))
//                          = log(|Sigma|^{-1/2}) = -sum(log(L[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// nx[1]        Dimension of the normal distribution
//
void
ldMVN1(double* log_dens,        double* work,
       const double* x,         const double* mu,  const double* Li,       
       const double* log_dets,  const int* nx);

void
ldMVN2(double* log_dens,        double* work,
       const double* x,         const double* mu,  const double* L,       
       const double* log_dets,  const int* nx);


/***** ***************************************************************************************** *****/
/***** Dist::ldMVN3                                                                              *****/
/***** ***************************************************************************************** *****/
//
// Log-density of the multivariate normal distribution N(mu, scale^{-1}*Q^{-1})
//
// log_dens[1]  Value of the log-density evaluated at x
// 
// work[nx]     Working array
//
// x[nx]        Point where the density is to be evaluated or generated random number
//
// mu[nx]       Mean of the normal distribution
//              * not needed when mu_nonZERO = 0 in which case it can be set to NULL
//
// Li[LT(nx)]   Cholesky decomposition of the precision matrix Q (lower triangle only), i.e., Q = Li %*% t(Li)
//              * array of length nx*(nx+1)/2
//
// log_dets[2]  log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j]))
//                          = log(|Sigma|^{-1/2}) = -sum(log(L[j,j]))
//              log_dets[1] = -nx*log(sqrt(2pi)) 
//
// sqrt_scale[1]        Scale factor for the covariance matrix
//
// log_sqrt_scale[1]    Log(sqrt(scale))
//
// nx[1]                Dimension of the normal distribution
//
void
ldMVN3(double* log_dens,          double* work,
       const double* x,           const double* mu,           
       const double* Li,          const double* log_dets,
       const double* sqrt_scale,  const double* log_sqrt_scale, 
       const int* nx);


/***** ************************************************************************************************* *****/
/***** Dist::dMVN1_R                                                                                     *****/
/***** Dist::rMVN1_R                                                                                     *****/
/***** ************************************************************************************************* *****/
//
// R wrappers to dMVN1 and rMVN1.
//
// log_dens[npoints]  Values of the (log-)density evaluated at x
// 
// Q[LT(nx)]          Precision matrix of the normal distribution (lower triangle only)
//                    on EXIT: re-written by its Cholesky decomposition
//                    * array of length nx*(nx+1)/2
//
// work[nx]           Working array
//
// err[1]             Error flag
//                    0 on EXIT if everything OK
//
// x[npoints*nx]      Points where the density is to be evaluated or generated random numbers
//                    * matrix npoints x nx in column major order where 1 column = 1 point
//
// unlog[1]           If not zero, then the value of the density is computed and not log-density
//
// mu[nx]             Mean of the normal distribution
//                    * not needed when mu_nonZERO = 0 in which case it can be set to NULL
//
// nx[1]              Dimension of the normal distribution
//
// mu_nonZERO[1]      If <> 0 then the argument 'mu' gives the mean of the normal distribution,
//                    otherwise it is assumed that the mean of the normal distribution is zero
//
// npoints[1]         Number of points where the (log-)density is to be evaluated
//                    or number of random numbers to generate
// 
#ifdef __cplusplus
extern "C" {
#endif

void
dMVN1_R(double* log_dens,  double* Q,              double* work,        int* err,
        const double* x,   const int* unlog,       const double* mu,
        const int* nx,     const int* mu_nonZERO,  const int* npoints);

void
rMVN1_R(double* x,         double* log_dens,  double* Q,              int* err,
        const double*mu,   const int* nx,     const int* mu_nonZERO,  const int* npoints);

#ifdef __cplusplus
}
#endif


/***** ************************************************************************************************* *****/
/***** Dist::rMVN2_R                                                                                     *****/
/***** ************************************************************************************************* *****/
//
// R wrappers to rMVN2.
//
// x[npoints*nx]      Generated random numbers
//                    * matrix npoints x nx in column major order where 1 column = 1 point
//
// mu[nx]             INPUT:   Canonical mean of the normal distribution
//                    OUTPUT:  Computed mean of the normal distribution
//
// log_dens[npoints]  Values of the (log-)density evaluated at x
// 
// Q[LT(nx)]          Precision matrix of the normal distribution (lower triangle only)
//                    on EXIT: re-written by its Cholesky decomposition
//                    * array of length nx*(nx+1)/2
//
// work[nx]           Working array
//
// err[1]             Error flag
//                    0 on EXIT if everything OK
//
// unlog[1]           If not zero, then the value of the density is computed and not log-density
//
// b[nx]              Canonical mean of the normal distribution
//
// nx[1]              Dimension of the normal distribution
//
// npoints[1]         Number of random numbers to generate
//                    
#ifdef __cplusplus
extern "C" {
#endif

void
rMVN2_R(double* x,      double* mu,         double* log_dens,       
        double* Q,      int* err,
        const int* nx,  const int* npoints);

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Dist ***/

#endif
