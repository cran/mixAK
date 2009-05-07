//
//  PURPOSE:   Truncated univariate normal distribution
//             * density, random numbers, ...
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/11/2007
//
//  FUNCTIONS:
//     * rTNorm1  13/11/2007:      Univariate truncated normal distribution, random numbers via inverse c.d.f. sampling
//
//     * rTNorm1_R  13/11/2007:    R wrapper for rTnorm1
//                    
// ==============================================================================================================================
//
#ifndef _DIST_TRUNCATED_NORM_H_
#define _DIST_TRUNCATED_NORM_H_

#include <R.h>
#include <Rmath.h>

namespace Dist{

const double N_limit = 8;         /** w = Phi^{-1}(p) cannot be computed if |w| > N_limit   **/
const double N_prob0 = 1e-15;     /** in R qnorm(1e-15)   = -7.94 (and is still not -Infty) **/
const double N_prob1 = 1-1e-15;   /** in R qnorm(1-1e-15) = 7.94 (and is still not Infty)   **/

/***** ***************************************************************************************** *****/
/***** Dist::rTNorm1                                                                             *****/
/***** ***************************************************************************************** *****/
//
// Inverse c.d.f. sampling
// X ~ TN(0, 1, a, b), then X = Phi^{-1}(U), where U ~ Unif(Phi(a), Phi(b))
// X ~ TN(mu, sigma^2, a, b), then X = Phi_{mu, sigma^2}^{-1}, where U ~ Unif(Phi_{mu, sigma^2}(a), Phi_{mu, sigma^2}(b))
//
// !!! Numerically, w = Phi^{-1}(p) cannot be computed if |w| > 8 !!!
//     In R: qnorm(pnorm(8))    = 7.991574
//           qnorm(pnorm(8.29)) = 8.209536
//           qnorm(pnorm(8.3))  = Infty
//
// That is why, if (after standardization), truncation interval is (>8, something), this function returns values equal to mu + sigma*8
//              if (after standardization), truncation interval is (something, <-8), this function returns values equal to mu - sigma*8
// So that in these extreme cases, we do not produce a sample from the truncated normal distribution!
//
// x[1]      OUTPUT: sampled value
//
// mu[1]     Mean of the untruncated normal distribution
//
// sigma[1]  Standard deviation of the untruncated normal distribution
//
// a[1]      Truncation limit 1
//           * used when trunc = 0, 2, 3
//
// b[1]      Truncation limit 2
//           * used only when trunc = 3
//
// trunc[1]  Type of truncation
//           * 0 = right-truncation, trunc. interval = (a, Infty) (b on INPUT is ignored)
//           * 1 = degenerated distribution P(X = a) = 1
//           * 2 = left-truncation, trunc. interval = (-Infty, a) (b on INPUT is ignored)
//           * 3 = interval-truncation, trunc. interval = (a, b)
//             it should hold: a < b (it is not checked in the C++ code!!!)
//           * 4 = no truncation, trunc. interval = (-Infty, Infty) (a and b on INPUT are ignored)
//
void
rTNorm1(double* x,  const double* mu,  const double* sigma,  const double* a,  const double* b,  const int* trunc);


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rTNorm1_R                                                                           *****/
/***** ***************************************************************************************** *****/
//
// x[nx]           OUTPUT: sampled value
//
// mu[1 or nx]     Mean(s) of the untruncated normal distribution
//
// sigma[1 or nx]  Standard deviation(s) of the untruncated normal distribution
//
// a[1 or nx]      Truncation limit(s) 1
//                 * used when trunc = 0, 2, 3
//
// b[1 or nx]      Truncation limit(s) 2
//                 * used only when trunc = 3
//
// trunc[1 or nx]  Type(s) of truncation
//                 * 0 = right-truncation, trunc. interval = (a, Infty) (b on INPUT is ignored)
//                 * 1 = degenerated normal distribution N(a, 0), i.e., P(X=a)=1
//                 * 2 = left-truncation, trunc. interval = (-Infty, a) (b on INPUT is ignored)
//                 * 3 = interval-truncation, trunc. interval = (a, b)
//                   it should hold: a < b (it is not checked in the C++ code!!!)
//                 * 4 = no truncation, trunc. interval = (-Infty, Infty) (a and b on INPUT are ignored)
//
// nx[1]           Number of points to sample
//
// mu_sigma_common[1]   If <> 0 then it is assumed that mu and sigma are the same for all sampled x's
//                      
// a_b_trunc_common[1]  If <> 0 then it is assumed that a, b and trunc are the same for all sampled x's
//
void
rTNorm1_R(double* x,      const double* mu,            const double* sigma,          const double* a,  const double* b,  const int* trunc,  
          const int* nx,  const int* mu_sigma_common,  const int* a_b_trunc_common);

#ifdef __cplusplus
}
#endif


}  /*** end of namespace Dist ***/

#endif
