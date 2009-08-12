//
//  PURPOSE:   (Multivariate) GLMM, update of the variances of continuous responses
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/07/2009
//
//  FUNCTIONS:  
//     * updateVars_eps  10/07/2009:  
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_VARS_EPS_H_
#define _GLMM_UPDATE_VARS_EPS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateVars_eps                                                                      *****/
/***** ***************************************************************************************** *****/
//
//  sigma[R_c]:     INPUT: whatsever
//                 OUTPUT: updated residual standard deviations
// 
//  Y[]:            values of continuous response (sorted as argument Y_c of GLMM_MCMC function)
//
//  eta[]:          values of linear predictors of continuous responses corresponding to Y
//
//  R[1]:           number of continuous responses
//
//  I[1]:           number of clusters
// 
//  n[R*I]:         numbers of observations for response and each cluster (sorted as argument n of GLMM_MCMC function)
//
//  N_s[R]:         total number of observations for each response (sum(n) = sum(N))
//
//  zeta[R]:        degrees of freedom of Wishart prior
//
//  gammaInv[R]:    inverted scale parameters of Wishart prior
//
void
updateVars_eps(double* sigma,  
               const double* Y,     const double* eta,  
               const int* R,        const int* I,            const int* n,  const int* N_s,
               const double* zeta,  const double* gammaInv);

}  /** end of namespace GLMM **/

#endif
