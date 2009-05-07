//
//  PURPOSE:   Normal mixture model, computation of the penalized expected deviance
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
//  FUNCTIONS:  
//     * NMix_PED  06/11/2008:  
//
// ====================================================================================================
//
#ifndef _NMIX_PENALIZED_EXPECTED_DEVIANCE_H_
#define _NMIX_PENALIZED_EXPECTED_DEVIANCE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"
#include "NMix_Utils.h"
#include "NMix_PED_coreUni.h"
#include "NMix_PED_coreMulti.h"
#include "Dist_TmixNorm.h"
#include "Dist_TmixMVN.h"

//namespace NMix{

#ifdef __cplusplus
extern "C" {
#endif

/***** ************************************************************************** *****/
/***** NMix_PED                                                                   *****/
/***** ************************************************************************** *****/
//
// PED[5]            INPUT:   whatsever
//                   OUTPUT:  PED[0] = Dbar = observed (expected) deviance
//                            PED[1] = p(opt) = total optimism (see p. 526 of Plummer, 2008)
//                                     obtained without reweighting the sample
//                            PED[2] = PED = PED[0] + PED[1]
//                            PED[3] = p(opt) obtained using the importance sampling
//                            PED[4] = PED[0] + PED[3]
//
// pm_indDevObs[n]   INPUT:   whatsever
//                   OUTPUT:  MCMC estimate of observed deviance for each observation,
//           i.e., pm_indDevObs[i] = average of 
//           0.5 * (-2 * log(sum_{j=1}^K(chain 1) w_j(chain 1) * phi(y_i | mu_j(chain 1), Sigma_j(chain 1))) 
//                  -2 * log(sum_{j=1}^K(chain 2) w_j(chain 2) * phi(y_i | mu_j(chain 2), Sigma_j(chain 2))))
//           over iterations
//
// pm_indpopt[n]      INPUT:   whatsever
//                   OUTPUT:  MCMC estimate of individual contributions to optimism obtained
//           without reweighting the sample
//
// pm_windpopt[n]     INPUT:   whatsever
//                   OUTPUT:  MCMC estimate of individual contributions to optimism obtained
//           using the importance sampling, see Plummer (2008)
//
// invalid_indDevObs[n]  INPUT:  whatsever
//                      OUTPUT:  number of iterations (in both chains) where log(f(y1))=Inf or log(f(y2))=Inf
//                               -> number between 0 and 2*M
//                               
// invalid_indpopt[n]    INPUT:  whatsever
//                      OUTPUT:  number of iterations where indDevObs = Inf or J(theta1, theta2) = +-Inf
//
// invalid_windpopt[n]   INPUT:  whatsever
//                      OUTPUT:  number of iterations where indDevObs = Inf or J(theta1, theta2) = +-Inf
//                               or weight for importance sampling = Inf
//
// sum_ISweight[n]   INPUT:   whatsever
//                   OUTPUT:  sum of relative weights for importance sampling for each observation
//           (see p. 530 of Plummer, 2008)
//
// y0[p, n]        data (matrix stored in column major order)
//                 * exactly observed data
//                 * right-censored observations                              
//                 * left-censored observations
//                 * lower limits of interval-censored observations
//                  
// y1[p, n]        data (matrix stored in column major order)
//                 * upper limits of interval-censored observations
//
// censor[p, n]    censoring indicator for each observation (matrix stored in column major order)
//                 * 0 = right-censored
//                 * 1 = observed
//                 * 2 = left-censored
//                 * 3 = interval-censored
//
// dimy[2]         dimension of the data
//                 * dimy[0]  = p = dimension of the response
//                 * dimy[+1] = n = number of observations
//
//
// chK1[1 or M]          Sampled numbers of mixture components (chain 1)
//                       chK1[0] = K if Krandom == 0
//
// chw1[sum(K1)]         Sampled weights (chain1)
//
// chmu1[p*sum(K1)]      Sampled means (chain1)
//
// chLi1[LT(p)*sum(K1)]  Cholesky factors of sampled inverse variances (chain1)
//
// chK2[1 or M]          Sampled numbers of mixture components (chain 2)
//                       chK2[0] = K if Krandom == 0
//
// chw2[sum(K2)]         Sampled weights (chain2)
//
// chmu2[p*sum(K2)]      Sampled means (chain2)
//
// chLi2[LT(p)*sum(K2)]  Cholesky factors of sampled inverse variances (chain2)
//
// M[1]       Lengths of the chains
//
// Kmax[1]    Maximal K
//
// Krandom[1]  If = 0 then it is assumed that K is constant
//
void
NMix_PED(double* PED,
         double* pm_indDevObs,    double* pm_indpopt,    double* pm_windpopt,
         int* invalid_indDevObs,  int* invalid_indpopt,  int* invalid_windpopt,
         double* sum_ISweight,    // double* ch_ISweight,
         int* err,
         const double* y0,         const double* y1,     const int* censor,    const int* dimy,
         const int* chK1,          const double* chw1,   const double* chmu1,  const double* chLi1,
         const int* chK2,          const double* chw2,   const double* chmu2,  const double* chLi2,
         const int* M,             const int* Kmax,      const int* Krandom,
         const double* Dens_ZERO,  const double* EMin);

#ifdef __cplusplus
}
#endif

//}  /*** end of namespace NMix ***/

#endif

