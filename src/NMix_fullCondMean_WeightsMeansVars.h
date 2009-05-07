//
//  PURPOSE:   Normal mixture model, mean of the full conditional distribution of mixture means and (inverse) variances
//             * needed to compute quantities for DIC_4 from Celeux, Forbes, Robert, Titterington (2006)
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/02/2008
//
//  FUNCTIONS:  * fullCondMean_MeansVars_NC  08/02/2008
//
// ====================================================================================================
//
#ifndef _NMIX_FULL_CONDITIONAL_MEAN_MEANS_VARS_H_
#define _NMIX_FULL_CONDITIONAL_MEAN_MEANS_VARS_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::fullCondMean_MeansVars_NC                                                           *****/
/***** ***************************************************************************************** *****/
//
// Compute means for the full conditional distribution of mixture means and inverse variances
// under the natural conjugate prior
//
// * code partially taken from NMix::updateMeansVars_NC
//
// fcm_weight[K]      Full conditional means for mixture weights
//
// logfcm_weight[K]   Logarithm of above
//
// fcm_mu[p, K]       Full conditional means for mixture means - E[mu_j|...]
// 
// ifcm_Q[LTp, K]     Inverted full conditional means for mixture inverse variances - (E[Q|...])^{-1}
//
// ifcm_L[LTp, K]     Cholesky decompositions of ifcm_Q's
//
// fcm_log_dets[2,K]  Log_dets based on normals where Sigma = (E[Q|...])^{-1}
//
// dwork[p]           Working space
//
// mixSumy[p, K]    mixSumy[,j] = sum_{i: r_i=j} y_i = n_j * ybar_j
//                  -> vector of 0 if n_j = 0
//
// mixBary[p, K]    mixBary[,j] = (1/n_j) * sum_{i: r_i=j} y_i = ybar_j
//                  -> vector of 0 if n_j = 0
//
// mixSS[LTp, K]    mixSS[,j] = sum_{i: r_i=j} (y_i - ybar_j) %*% t(y_i - ybar_j)
//                  -> matrix of 0 if n_j = 0
//
// p[1]             dimension of the normal distribution
//
// n[1]             number of observations
//
// K[1]             current number of components
//
// Q[NULL]          nowhere used
//                  * it is here for prototype compatibility with NMix::updateMeansVars_IC
//
// delta[1]         parameter of the Dirichlet prior on weights
//
// c[K]             prior precisions c_1, ..., c_K
//
// xi[p,K]          prior means 
//
// c_xi[p,K]        prior means scaled by prior precisions
//                  * c_xi[,j] = c_j * xi_j
//
// Dinv[NULL]       nowhere used
//                  * it is here for prototype compatibility with NMix::updateMeansVars_IC
//
// Dinv_xi[NULL]    nowhere used
//                  * it is here for prototype compatibility with NMix::updateMeansVars_IC
//
// zeta[1]          prior degrees of freedom of the Wishart distribution
//
// XiInv[LT(p)]     inverse of the prior scale matrix of the Wishart distribution
//                  * a priori E(Sigma_j^{-1}) = zeta * Xi
//
void
fullCondMean_WeightsMeansVars_NC(double* fcm_weight,     double* logfcm_weight,  double* fcm_mu,         
                                 double* ifcm_Q,         double* ifcm_L,         double* fcm_log_dets,  
                                 double* dwork,          int* err,
                                 const double* mixSumy,  const double* mixBary,  const double* mixSS,
                                 const int* mixN,        const int* p,           const int* n,         const int* K,         const double* Q,
                                 const double* delta,    const double* c,        const double* xi,     const double* c_xi,   
                                 const double* Dinv,     const double* Dinv_xi,  const double* zeta,   const double* XiInv);

}  /*** end of namespace NMix ***/

#endif
