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
//             * random numbers, ...
//
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111208  created
//
//  FUNCTIONS:
//     * rMVT  20111208:   Random numbers from MVT distribution
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


#ifdef __cplusplus
extern "C" {
#endif

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
void
rMVT1_R(double* x,         
        double* log_dens,  
        double* Q,              
        int*    err,
        const double* nu,
        const double* mu,  
        const int*    nx,
        const int*    npoints);

#ifdef __cplusplus
}
#endif

}  // end of namespace Dist

#endif
