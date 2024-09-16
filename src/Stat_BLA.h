//
//  PURPOSE:   Best linear approximation (theoretical least squares)
//             * useful for computation of conditional normal distributions
//
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/11/2007
//
//  FUNCTIONS:
//     * BLA  19/11/2007:      Compute regression coefficients for:  
//                             x[0]  = beta[0,0] + beta[1,0]*x[1] + ... + beta[p-1,0]*x[p-1]
//                             x[1]  = beta[0,1] + beta[1,1]*x[0] + ... + beta[p-1,1]*x[p-1]
//                                      ...
//                             x[p-1]= beta[0,p-1] + beta[1,p-1]*x[0] + ... + beta[p-1,p-1]*x[p-2]
//
//                    
// ==============================================================================================================================
//
#ifndef _STAT_BEST_LINEAR_APPROXIMATION_H_
#define _STAT_BEST_LINEAR_APPROXIMATION_H_

#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace Stat{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Stat::BLA                                                                                 *****/
/***** ***************************************************************************************** *****/
//
// Compute regression coefficients for:  x[0]  = beta[0,0] + beta[1,0]*x[1] + ... + beta[p-1,0]*x[p-1]
//                                       x[1]  = beta[0,1] + beta[1,1]*x[0] + ... + beta[p-1,1]*x[p-1]
//                                       ...
//                                       x[p-1]= beta[0,p-1] + beta[1,p-1]*x[0] + ... + beta[p-1,p-1]*x[p-2]
//
// beta[p x p]        Computed regression coefficients (in columns)
//                    * beta[1:(p-1), i] = Sigma[-i,-i]^{-1} %*% Sigma[i,-i]
//                    * beta[0, i]       = mu[i] - t(beta[1:(p-1), i]) %*% mu[-i]
//
// sigmaR2[p]         Residual variances for each of p regressions
//                    sigmaR2[i] = Sigma[i,i] - t(beta[1:(p-1), i]) %*% Sigma[-i,-i] %*% beta[1:(p-1), i]
//                               = Sigma[i,i] - Sigma[i,-i] %*% Sigma[-i,-i]^{-1} %*% Sigma[-i,i]
//
// L[LT(p-1)]         Working space to store lower triangles of matrices L_i = Cholesky decomposition of Sigma[-i,-i] (i=0, ..., p-1)
//
// err[1]             Error flag, equal to 0 if OK
//
// mu[p]              Vector of means
//
// Sigma[LT(p)]       Lower triangle of var(X)
//
// p[1]               Dimension of the random vector (it should be >= 2 -> NOT CHECKED)
//
void
BLA(double* beta,      
    double* sigmaR2,      
    double* L,     
    int*    err,  
    const double* mu,  
    const double* Sigma,  
    const int*    p);

#ifdef __cplusplus
}
#endif

}    /*** end of namespace Stat ***/

#endif
