//
//  PURPOSE:   Truncated multivariate normal distribution
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/11/2007
//
//  FUNCTIONS:
//     * rTMVN1  20/11/2007:      Multivariate truncated normal distribution, random numbers via Gibbs sampling
//                                as described in Geweke (1991)
//
//     * rTMVN2  19/11/2007:      Multivariate truncated normal distribution, random numbers via Gibbs sampling
//                                as described in Rodriguez-Yam, Davis, Scharf (submitted in 2004)
//                                - NOT FINISHED
//
//     * rTMVN1_R  20/11/2007:    R wrapper for rTMVN1
//                    
// ==============================================================================================================================
//
#ifndef _DIST_TRUNCATED_MVN_H_
#define _DIST_TRUNCATED_MVN_H_

#include <R.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"
#include "Dist_TNorm.h"
#include "Stat_BLA.h"

namespace Dist{


/***** ***************************************************************************************** *****/
/***** Dist::rTMVN1                                                                              *****/
/***** ***************************************************************************************** *****/
//  
// Perform one MCMC scan of the Gibbs algorithm described in Geweke (1991).
// The stationary distribution of so created Markov chain is the truncated (multivariate) normal distribution.
//
// We wish to sample from the TMVN(mu, Sigma, -c0 <= x <= c1), where
//   * mu is the mean vector of length p
//
//   * Sigma is pxp positive definite covariance matrix
//
//   * c0 is the vector of length r=p, its components are allowed to be equal to Infty
//
//   * c1 is the vector of length r=p, its components are allowed to be equal to Infty
//
//
// x[p]    INPUT:   Initial value to start MCMC or the previous state of the Markov chain
//         OUTPUT:  Updated value of the Markov chain
// 
// beta[p x p]        Regression coefficients for regressions x[i] ~ x[-i] (in columns)
//                    * beta[1:(p-1), i] = Sigma[-i,-i]^{-1} %*% Sigma[i,-i]
//                    * beta[0, i]       = mu[i] - t(beta[1:(p-1), i]) %*% mu[-i]
//
// sigmaR2[p]         Residual variances for each of p regressions
//                    sigmaR2[i] = Sigma[i,i] - t(beta[1:(p-1), i]) %*% Sigma[-i,-i] %*% beta[1:(p-1), i]
//                               = Sigma[i,i] - Sigma[i,-i] %*% Sigma[-i,-i]^{-1} %*% Sigma[-i,i]
//
// a[p]         Truncation limits 1, see trunc below for more details
//
// b[p]         Truncation limits 2, see trunc below for more details
//
// trunc[p]     Types of truncation
//              * 0:  a_i <= x_i < Infty,     i.e. c_i = -a_i,  c_{p+i} = Infty
//              * 1:  a_i <= x_i <= a_i,      i.e. c_i = -a_i,  c_{p+i} = a_i
//              * 2:  -Infty < x_i <= a_i,    i.e. c_i = Infty, c_{p+i} = a_i
//              * 3:  a_i <= x_i <= b_i,      i.e. c_i = -a_i,  c_{p+i} = b_i
//              * 4:  -Infty <= x_i <= Infty, i.e. c_i = Infty, c_{p+i} = Infty
//
// p[1]         Dimension of the normal distribution
//
void
rTMVN1(double* x,  
       const double* beta,  const double* sigmaR2,  
       const double* a,     const double* b,        const int* trunc,  const int* p);


/***** ***************************************************************************************** *****/
/***** Dist::rTMVN2                                                                              *****/
/***** ***************************************************************************************** *****/
//  
// Perform one MCMC scan of the Gibbs algorithm described in Rodriguez-Yam, Davis, Scharf (submitted in 2004).
// The stationary distribution of so created Markov chain is the truncated (multivariate) normal distribution.
//
// REMARK:  It is still unclear to me what to do when some inequalities are in fact equalities.
//          Also, I am afraid that the algorithm is not too good if there are some small truncation intervals.
//
// We wish to sample from the TMVN(mu, Sigma, B*x <= c), where
//   * mu is the mean vector of length p
//
//   * Sigma is pxp positive definite covariance matrix
//
//   * B is in general rxp matrix whose rows are allowed to be linearly dependent
//     This version assumes that t(B) = (-I_p I_p), that is the constraints are
//     x_1 >= -c_1, ..., x_p >= -c_p, 
//     x_1 <= c_{p+1}, ..., x_p <= c_{p+p} 
//
//   * c is the vector of length r=2*p, its components are allowed to be equal to Infty
//
//   * Let Q = Sigma^{-1}
//
//   * Let Q = Li %*% t(Li) is Cholesky decomposition of Q
//
//   * Sigma = L %*% t(L) is Cholesky decomposition of Sigma
//
//   * Let A be pxp matrix such that A %*% Sigma %*% t(A) = I_p
//     NOTE: You can take A = L^{-1} (lower triangular matrix) or A = t(Li) (upper triangular matrix)
//           It is assumed in this version that A is either L^{-1} or t(Li)
//
//   * Let G = A^{-1}
//
//
// x[p]    INPUT:   Initial value to start MCMC or the previous state of the Markov chain
//         OUTPUT:  Updated value of the Markov chain
// 
// dwork[p + p] Working array
//
// alpha[p]     Transformed mean of the untruncated normal distribution, i.e., alpha = A %*% mu
//
// Ainv[LT(p)]  Either L or Li^{-1} (see whichA for more details) stored in a packed form
//
// whichA[1]    Indicates whether A = L^{-1} or A = t(Li)
//              * whichA <> 0 => A = t(Li) and array Ainv contains lower triangle of Li^{-1}
//              * whichA == 0 => A = L^{-1} and array Ainv contains lower triangle of L
//
// a[p]         Truncation limits 1, see trunc below for more details
//
// b[p]         Truncation limits 2, see trunc below for more details
//
// trunc[p]     Types of truncation
//              * 0:  a_i <= x_i < Infty,     i.e. c_i = -a_i,  c_{p+i} = Infty
//              * 1:  a_i <= x_i <= a_i,      i.e. c_i = -a_i,  c_{p+i} = a_i
//              * 2:  -Infty < x_i <= a_i,    i.e. c_i = Infty, c_{p+i} = a_i
//              * 3:  a_i <= x_i <= b_i,      i.e. c_i = -a_i,  c_{p+i} = b_i
//              * 4:  -Infty <= x_i <= Infty, i.e. c_i = Infty, c_{p+i} = Infty
//
// p[1]         Dimension of the normal distribution
//
void
rTMVN2(double* x,             double* dwork,
       const double* alpha,   const double* Ainv,   const int*  whichA,
       const double* a,       const double* b,      const int* trunc,    const int* p);


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rTMVN1_R                                                                            *****/
/***** ***************************************************************************************** *****/
//
// R wrapper to rTMVN1_R
//
// x[npoints*p]    Generated Markov chain
//                 * matrix npoints x p in column major order where 1 column = 1 point
//
// beta[p*p]       Regression coefficients computed using Dist::BLA
//
// sigmaR2[p]      Residual variances computed using Dist::BLA
//
// L[LT(p-1)]      Working space for Dist::BLA
//
// err[1]          Error flag
//                 0 on EXIT if everything OK
//
// mu[p]           Mean of the normal distribution before truncation
//
// Sigma[LT(p)]    Lower triangle of the covariance matrix of the normal distribution before truncation
//
// a[p]            Truncation limits 1 (see trunc below)
//
// b[p]            Truncation limits 2 (see trunc below)
//
// trunc[p]        Types of truncation in each margin (see rTMVN1)
//
// p[1]            Dimension of the normal distribution
//
// npoints[1]      Number of random numbers to generate
//
void
rTMVN1_R(double* x,            double* beta,         double* sigmaR2,    
         double* L,            int* err,
         const double* xinit,  const double* mu,     const double* Sigma,          
         const double* a,      const double* b,      const int* trunc,  
         const int* p,         const int* npoints);

#ifdef __cplusplus
}
#endif


}  /*** end of namespace Dist ***/

#endif
