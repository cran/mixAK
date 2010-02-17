//
//  PURPOSE:   Normal mixture model, 
//             computation of 
//             Pr_y[j, i] = w_j * phi(y_i | mu_j, Sigma_j), i=0,...,n-1, j=0,...,K-1
//             and related cumulative sums
//             cum_Pr_y[j, i] = sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l), i=0,...,n-1, j=0,...,K-1
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/02/2010
//
//  FUNCTIONS:  
//     * Pr_y_and_cum_Pr_y  10/02/2010:  
//
// ===================================================================================
//
#ifndef _NMIX_PR_Y_AND_CUM_PR_Y_H_
#define _NMIX_PR_Y_AND_CUM_PR_Y_H_

#include <R.h>

#include "AK_Basic.h"
#include "Dist_MVN.h"

namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** NMix::Pr_y_and_cum_Pr_y                                                                   *****/
/***** ***************************************************************************************** *****/
//
//  REMARK:  Pr_y     are re-scaled to sum-up to 1
//           cum_Pr_y are cumulative sums which are NOT re-scaled to give a total of 1
//
//  Pr_y[K, n]      INPUT:  whatsever  
//                  OUTPUT: Pr_y[j, i] = w_j * phi(y_i | mu_j, Sigma_j), i=0,...,n-1, j=0,...,K-1
//
//  cum_Pr_y[K, n]  INPUT:  whatsever  
//                  OUTPUT: cum_Pr_y[j, i] \propto sum_{l=1}^j w_l * phi(y_i | mu_l, Sigma_l), i=0,...,n-1, j=0,...,K-1
//
//  dwork[p]        working array for dMVN
//
//  y[p, n]         data, y[, i] = y_i
//  
//  p[1]            dimension of the data
//
//  n[1]            sample size
//
//  logw[K]         logarithm of mixture weights
//
//  mu[p, K]        mixture means
//
//  Li[LT(p), K]    Cholesky decompositions of mixture inverse variances
//
//  log_dets[2, K]  input for dMVN
//
//  K[1]            number of mixture components
//
void
Pr_y_and_cum_Pr_y(double* Pr_y,
                  double* cum_Pr_y,
                  double* dwork,
                  const double* y,
                  const int*    p,
                  const int*    n,
                  const double* logw,
                  const double* mu,
                  const double* Li,
                  const double* log_dets,
                  const int*    K);

#ifdef __cplusplus
}
#endif

}  // end of namespace NMix

#endif

