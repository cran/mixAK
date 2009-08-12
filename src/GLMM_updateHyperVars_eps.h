//
//  PURPOSE:   (Multivariate) GLMM, update of the hyperparameters of the variances
//             of continuous responses
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/07/2009
//
//  FUNCTIONS:  
//     * updateHyperVars_eps  10/07/2009:  
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_HYPER_VARS_EPS_H_
#define _GLMM_UPDATE_HYPER_VARS_EPS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateHyperVars_eps                                                                 *****/
/***** ***************************************************************************************** *****/
//
//
//  ASSUMPTION:  For each response variable:
//               sigma^{-1} ~ Wishart(zeta, gamma) = Gamma(zeta/2, gamma^{-1}/2)
//               gamma^{-1} ~ Gamma(g, h)
//
//  ---------------------------------------------------------------------------------------------------
//
//  gammaInv[R]      INPUT:  whatsever   
//                  OUTPUT:  inverse gamma's
//
//  sigma[R]         residual standard deviations
//  R[1]             number of std. deviations
//  zeta[R]          prior degrees of freedom of the Wishart distribution
//  g[R]             parameter of gamma prior
//  h[R]             parameter of gamma prior
//
void
updateHyperVars_eps(double* gammaInv,  
                    const double* sigma,  const int* R,
                    const double* zeta,   const double* g,  const double* h);

}  /** end of namespace GLMM **/

#endif
