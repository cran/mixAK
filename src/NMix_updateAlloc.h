//
//  PURPOSE:   Normal mixture model, update of the allocation variables
//             and related quantities
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2007
//
//  FUNCTIONS:  
//     * updateAlloc  06/11/2007:  Update of the allocations
//                    21/12/2007:  Validated in R, some crucial bugs fixed
//                    12/02/2008:  argument cum_Pr and related added
//
// ===================================================================================
//
#ifndef _NMIX_UPDATE_ALLOC_H_
#define _NMIX_UPDATE_ALLOC_H_

#include <R.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"
#include "Dist_Discrete.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::updateAlloc                                                                         *****/
/***** ***************************************************************************************** *****/
//
// r[dimy[1]]           INPUT:   whatsever
//                      OUTPUT:  updated allocations
//
// mixN[Kmax]           INPUT:   whatsever
//                      OUTPUT:  updated numbers of observations belonging to each component
//                               * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// rInv[KMax][n]        INPUT:   whatsever
//                      OUTPUT:  indeces of columns of y indicating observations belonging to each component
//
// cum_Pr[K, n]         INPUT:   if cum_Pr_done = false then whatsever
//                               if cum_Pr_done = true then cum_Pr[j, i] = sum_{l=0}^j w_l * phi(y_i | mu_l, Sigma_l)
//                      OUTPUT:  always cum_Pr[j, i] = C * sum_{l=0}^j w_l * phi(y_i | mu_l, Sigma_l)
//
// dwork_ldMVN[p]       working array
//
// y[p, n]              observations (in columns)
//
// p[1]                 dimension of the response
//
// n[1]                 number of observations
//
// logw[K]              mixture log-weights
//
// mu[p, K]             mixture means (in "columns")
//                      * mu_j = mu[,j]
//
// Li[LT(p), K]         Cholesky decompositions (lower triangles only) of mixture precision matrices (in columns)
//                      * Sigma_j^{-1} = Li[,j] %*% t(Li[,j])
//
// log_dets[2, K]       factors to compute log-density of normal distributions
//                      * log_dets[0, j] = log(|Sigma[j]|^{-1/2}) = sum(log(Li_{j}[l,l]))
//                      * log_dets[1, j] = -p*log(sqrt(2*pi))
//
// K[1]                 current number of components
//
// cum_Pr_done[1]       true/false, see above
//
void 
updateAlloc(int* r,   
            int* mixN,   
            int** rInv,   
            double* cum_Pr,   
            double* dwork_ldMVN,
            const double* y,      
            const int* p,       
            const int* n,
            const double* logw,   
            const double* mu,   
            const double* Li,   
            const double* log_dets,  
            const int* K,  
            const bool* cum_Pr_done);

}   /** end of namespace NMix **/

#endif

