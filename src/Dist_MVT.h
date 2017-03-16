//
//  PURPOSE:   Multivariate Student t distribution MVT_p(mu, Sigma, nu)
//             defined as T = mu + L*U*sqrt(nu/V),
//             where mu in Real^p
//                   Sigma = L %*% t(L) is positive definite p x p matrix and L is lower triangular
//                   nu > 0 are degrees of freedom
//                   U ~ N_p(0, I)
//                   V ~ chisq_nu independent of U
//                   or in more general V ~ Gamma(nu/2, 1/2) independent of U
//                               
//             modus(T) = mu
//             E(T)     = mu for nu > 1
//             var(T)   = (nu/(nu - 2)) * Sigma for nu > 2
//
//             * random numbers, density
//
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111208  functions for random numbers generation written
//             20120124  functions for evaluation of a density added
//
//  FUNCTIONS:
//     * rMVT1           20111208:   Random numbers from MVT distribution  
//     * ldMVT1          20120124:   Density of the MVT distribution
//     * deriv_ldMVT_x   20120124
//     * rMVT1_R         20111208:   Wrapper of rMVT1 to R
//     * dMVT1_R         20120124:   Wrapper of the ldMVT1 to R
//
// 
// ===================================================================================================
//
#ifndef _DIST_MVT_H_
#define _DIST_MVT_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>       

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rMVT1                                                                               *****/
/***** ***************************************************************************************** *****/
//
//  x[nx]         INPUT:  whatsever
//               OUTPUT:  sampled value
//
//  log_dens[1]   INPUT:  whatsever
//               OUTPUT:  value of the log-density at sampled x
//
//  nu[1]        degrees of freedom of the MVT distribution
//
//  mu[nx]       location parameter of the MVT distribution
//
//  Li[LT(nx)]   lower triangle of lower triangular matrix Li, where Q = Sigma^{-1} = Li %*% t(Li)
//
//  log_dets[1]  log_dets[0] = log(|Q|^{1/2}) = log(|Sigma|^{-1/2}) = sum(log(Li[j,j]))
//               log_dets[1] = lgamma((nu + nx)/2) - lgamma(nu/2) - (nx/2)*(log(nu) + log(pi)) 
//
//  nx[1]        dimension of the MVT distribution 
//
void
rMVT1(double*       x,
      double*       log_dens,
      const double* nu,
      const double* mu,
      const double* Li,
      const double *log_dets,
      const int*    nx);


/***** ***************************************************************************************** *****/
/***** Dist::rMVT1                                                                               *****/
/***** ***************************************************************************************** *****/
//
//  log_dens[1]   INPUT:  whatsever
//               OUTPUT:  value of the log-density evaluated at x
//
//  work[nx]     working array
//
//  x[nx]        point to evaluate the density
//
//  nu[1]        degrees of freedom of the MVT distribution
//
//  mu[nx]       location parameter of the MVT distribution
//
//  Li[LT(nx)]   lower triangle of lower triangular matrix Li, where Q = Sigma^{-1} = Li %*% t(Li)
//
//  log_dets[1]  log_dets[0] = log(|Q|^{1/2}) = log(|Sigma|^{-1/2}) = sum(log(Li[j,j]))
//               log_dets[1] = lgamma((nu + nx)/2) - lgamma(nu/2) - (nx/2)*(log(nu) + log(pi)) 
//
//  nx[1]        dimension of the MVT distribution 
//
void
ldMVT1(double*       log_dens,
       double*       work,
       const double* x,    
       const double* nu,
       const double* mu,
       const double* Li,
       const double *log_dets,
       const int*    nx);


/***** ***************************************************************************************** *****/
/***** Dist::deriv_ldMVT_x                                                                       *****/
/***** ***************************************************************************************** *****/
//
//  Function to calculate the first and the second derivatives w.r.t. x
//  of the (-(nu + p)/2) * log(1 + (t(x-mu) %*% Sigma^{-1} %*% (x-mu))/nu)
//  part of the log-density of the MVT distribution.
//  Needed when evaluating the Laplace approximation of the log-likelihood in a model
//  with random effects having the MVT distribution.
//
//  gradient[nx]     INPUT:  whatsever
//                  OUTPUT:  calculated first derivatives
//
//  Hessian[LT(nx)]  INPUT:  whatsever
//                  OUTPUT:  calculated second derivatives (lower triangle)
//
//  x[nx]        point to evaluate the derivatives
//
//  nu[1]        degrees of freedom of the MVT distribution
//
//  mu[nx]       location parameter of the MVT distribution
//
//  Q[LT(nx)]    lower triangle of a symmetric matrix Q = Sigma^{-1}
//
//  Li[LT(nx)]   lower triangle of lower triangular matrix Li, where Q = Sigma^{-1} = Li %*% t(Li)
//
//  nx[1]        dimension
//
void
deriv_ldMVT_x(double*       gradient,
              double*       Hessian,
              const double* x,    
              const double* nu,
              const double* mu,
              const double* Q,
              const double* Li,
              const int*    nx);



/***** ***************************************************************************************** *****/
/***** Dist::rMVT1_R                                                                           *****/
/***** ***************************************************************************************** *****/
//
//  x[nx, npoints]    INPUT:  whatsever
//                   OUTPUT:  sampled value
//
//  log_dens[1]   INPUT:  whatsever
//               OUTPUT:  value of the log-density at sampled x
//
//  Q[LT(nx)]     INPUT:  precission matrix Q = Sigma^{-1/2} (its lower traingle)
//               OUTPUT:  Cholesky decomposition of Q
//
//  err[1]       error flag
//
//  nu[1]        degrees of freedom of the MVT distribution
//
//  mu[nx]       location parameter of the MVT distribution
//
//  nx[1]        dimension of the MVT distribution 
//
//  npoints[1]   number of points to sample from MVT
#ifdef __cplusplus
extern "C" {
#endif

void
rMVT1_R(double* x,         
        double* log_dens,  
        double* Q,              
        int*    err,
        const double* nu,
        const double* mu,  
        const int*    nx,
        const int*    npoints);


/***** ***************************************************************************************** *****/
/***** Dist::dMVT1_R                                                                           *****/
/***** ***************************************************************************************** *****/
//
//  log_dens[1]   INPUT:  whatsever
//               OUTPUT:  value of the log-density at x values
//
//  Q[LT(nx)]     INPUT:  precission matrix Q = Sigma^{-1/2} (its lower traingle)
//               OUTPUT:  Cholesky decomposition of Q
//
//  work[nx]     working array
//
//  err[1]       error flag
//
//  x[nx, npoints]   values to evaluate the density
//
//  unlog[1]     0 if log-density should be returned
//               1 if density should be returned
//
//  nu[1]        degrees of freedom of the MVT distribution
//
//  mu[nx]       location parameter of the MVT distribution
//
//  nx[1]        dimension of the MVT distribution 
//
//  npoints[1]   number of points to evaluate the density
  //
void
dMVT1_R(double*       log_dens,  
        double*       Q,              
        double*       work,       
        int*          err,
        const double* x,   
        const int*    unlog,
        const double* nu,
        const double* mu,      
        const int*    nx,     
        const int*    npoints);

#ifdef __cplusplus
}
#endif

}  // end of namespace Dist

#endif
