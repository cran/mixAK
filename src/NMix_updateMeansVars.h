//
//  PURPOSE:   Normal mixture model, update of the mixture means and variances
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/11/2007
//
//  FUNCTIONS:  
//     * updateMeansVars_NC  12/11/2007:   Update of the mixture means and variances when a natural conjugate prior
//                                         (Normal-Wishart) is assumed
//                           21/12/2007:   Partially validated in R, some crucial bugs fixed
//
//     * updateMeansVars_IC  13/02/2008:   Update of the mixture means and variances when an independent conjugate prior is assumed
//
// ====================================================================================================
//
#ifndef _NMIX_UPDATE_MEANS_VARS_H_
#define _NMIX_UPDATE_MEANS_VARS_H_

#include <R.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"
#include "Dist_Wishart.h"

#include "NMix_Utils.h"
#include "NMix_orderComp.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateMeansVars_NC                                                                  *****/
/***** ***************************************************************************************** *****/
//
// Natural conjugate prior
//
// Prior:             p(mu_j, Sigma_j^{-1}) = p(mu_j | Sigma_j) * p(Sigma_j^{-1})
//                    mu_j | Sigma_j ~ N(xi_j, c_j^{-1}*Sigma_j)
//                      Sigma_j^{-1} ~ W(zeta, Xi)
//                    
// Full conditionals:  mu_j | Sigma_j^{-1}, ... ~ N(m_j, (c_j + n_j)^{-1}*Sigma_j), where
//                            m_j      = (n_j + c_j)^{-1} * (mixSumy_j + c_j*xi_j)
//
//                     Sigma_j^{-1} | ... ~ W(zeta + n_j, (Xi^{-1} + S_j + ((c_j*n_j)/(c_j + n_j))*(yBar_j - xi_j)*t(yBar_j - xi_j))^{-1}),
//                            where yBar_j = n_j^{-1}*mixSumy_j
//                                  S_j    = sum_{i: r_i=j} (y_i - yBar_j)*t(y_i - yBar_j)
//
// mu[p, K]             INPUT:   whatsever
//                      OUTPUT:  updated mixture means (mu_j)
//
// Q[LT(p), K]          INPUT:   whatsever
//                      OUTPUT:  updated mixture inverse variances (Q_j = Sigma_j^{-1})
//
// Li[LT(p), K]         INPUT:   whatsever
//                      OUTPUT:  Cholesky decompositions (lower triangles only) of updated mixture precision matrices (in columns)
//                               * Q_j = Sigma_j^{-1} = Li[,j] %*% t(Li[,j])
//
// Sigma[LT(p), K]      INPUT:   whatsever
//                      OUTPUT:  mixture variances
//
// log_dets[2, K]       INPUT:   whatsever
//                      OUTPUT:  factors to compute log-density of the component normal distributions
//                               * log_dets[0, j] = log(|Sigma[j]|^{-1/2}) = sum(log(Li_{j}[l,l]))
//                               * log_dets[1, j] = -p*log(sqrt(2*pi))
//
// order[K]             INPUT:   whatsever
//                      OUTPUT:  order indeces of the mixture components after update
//
// rank[K]              INPUT:   whatsever
//                      OUTPUT:  rank indeces of the mixture components after update
//
// dwork[p + LT(p) + 2 + p + LT(p) + 2*p*p + K] working array
//
// err[1]               OUTPUT:  error flag
//                               (normally, there is no real reason for failure here...)
//
// y[p, n]
//
// r[n]
//
// mixN[K]
//
// p[1]
//
// n[1] 
//
// K[1]                 current number of components
//
// p[1]                 dimension of the normal distribution
//
// c[K]                 prior precisions c_1, ..., c_K
//
// xi[p,K]              prior means 
//
// c_xi[p,K]            prior means scaled by prior precisions
//                      * c_xi[,j] = c_j * xi_j
//
// Dinv[NULL]           nowhere used
//                      * it is here for prototype compatibility with NMix::updateMeansVars_IC
//
// Dinv_xi[NULL]        nowhere used
//                      * it is here for prototype compatibility with NMix::updateMeansVars_IC
//
// zeta[1]              prior degrees of freedom of the Wishart distribution
//
// XiInv[LT(p)]         inverse of the prior scale matrix of the Wishart distribution
//                      * a priori E(Sigma_j^{-1}) = zeta * Xi
//
void
updateMeansVars_NC(double* mu,          double* Q,              double* Li,          double* Sigma,
                   double* log_dets,    int* order,             int* rank,           double* dwork,       int* err,
                   const double* y,     const int* r,           const int* mixN,     const int* p,        const int* n,
                   const int* K,        const double* c,        const double* xi,    const double* c_xi,  
                   const double* Dinv,  const double* Dinv_xi,  const double* zeta,  const double* XiInv);


/***** ***************************************************************************************** *****/
/***** NMix::updateMeansVars_IC                                                                  *****/
/***** ***************************************************************************************** *****/
//
// Independent conjugate prior
//
// Prior:             p(mu_j, Sigma_j^{-1}) = p(mu_j) * p(Sigma_j^{-1})
//                              mu_j ~ N(xi_j, D_j)
//                      Sigma_j^{-1} ~ W(zeta, Xi)
//                    
// Full conditionals:  mu_j | ... ~ N(m_j, S_j), where
//                            S_j^{-1} = n_j*Sigma_j^{-1} + D_j^{-1}
//                            m_j      = S_j * (Sigma_j^{-1}*mixSumy_j + D_j^{-1}*xi_j)
//
//                     Sigma_j^{-1} | ... ~ W(zeta + n_j, (Xi^{-1} + sum_{i: r_i=j}(y_i - mu_j)*t(y_i - mu_j))^{-1}),
//                            where yBar_j = n_j^{-1}*mixSumy_j
//
// mu[p, K]             INPUT:   whatsever
//                      OUTPUT:  updated mixture means (mu_j)
//
// Q[LT(p), K]          INPUT:   whatsever
//                      OUTPUT:  updated mixture inverse variances (Q_j = Sigma_j^{-1})
//
// Li[LT(p), K]         INPUT:   whatsever
//                      OUTPUT:  Cholesky decompositions (lower triangles only) of updated mixture precision matrices (in columns)
//                               * Q_j = Sigma_j^{-1} = Li[,j] %*% t(Li[,j])
//
// Sigma[LT(p), K]      INPUT:   whatsever
//                      OUTPUT:  mixture variances
//
// log_dets[2, K]       INPUT:   whatsever
//                      OUTPUT:  factors to compute log-density of the component normal distributions
//                               * log_dets[0, j] = log(|Sigma[j]|^{-1/2}) = sum(log(Li_{j}[l,l]))
//                               * log_dets[1, j] = -p*log(sqrt(2*pi))
//
// order[K]             INPUT:   whatsever
//                      OUTPUT:  order indeces of the mixture components after update
//
// rank[K]              INPUT:   whatsever
//                      OUTPUT:  rank indeces of the mixture components after update
//
// dwork[p + LT(p) + 2 + p + LT(p) + 2*p*p + K] working array
//
// err[1]               OUTPUT:  error flag
//                               (normally, there is no real reason for failure here...)
//
// y[p, n]
//
// r[n]
//
// mixN[K]
//
// p[1]
//
// n[1] 
//
// K[1]                 current number of components
//
// p[1]                 dimension of the normal distribution
//
// c[K]                 prior precisions c_1, ..., c_K
//
// xi[NULL]             nowhere used
//                      * it is here for prototype compatibility with NMix::updateMeansVars_NC
//
// c_xi[NULL]           nowhere used
//                      * it is here for prototype compatibility with NMix::updateMeansVars_NC
//
// Dinv[LT(p),K]        prior inverse variances for mixture means
//
// Dinv_xi[p,K]         D_j^{-1} * xi_j
//
// zeta[1]              prior degrees of freedom of the Wishart distribution
//
// XiInv[LT(p)]         inverse of the prior scale matrix of the Wishart distribution
//                      * a priori E(Sigma_j^{-1}) = zeta * Xi
//
void
updateMeansVars_IC(double* mu,          double* Q,              double* Li,          double* Sigma,
                   double* log_dets,    int* order,             int* rank,           double* dwork,       int* err,
                   const double* y,     const int* r,           const int* mixN,     const int* p,        const int* n,
                   const int* K,        const double* c,        const double* xi,    const double* c_xi,  
                   const double* Dinv,  const double* Dinv_xi,  const double* zeta,  const double* XiInv);

}   /** end of namespace NMix **/

#endif

