//
//  PURPOSE:   Compute a mean and decomposition of the inverted covariance matrix of the
//             normal approximation based on one Newton-Raphson step
//             - computation is done using the LS solution and QR decomposition
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/04/2010
//
//  FUNCTIONS:  
//     *   11/04/2010:  MCMC:Moments_NormalApprox_QR (PROTOTYPE 1)
//     *   14/04/2010:  MCMC:Moments_NormalApprox_QR (PROTOTYPE 2)
//
// =================================================================================
//
#ifndef _MCMC_MOMENTS_NORMALAPPROX_QR_H_
#define _MCMC_MOMENTS_NORMALAPPROX_QR_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>                // contains declaration for dqrdc2 (QR decomposition) and dqrls (LS solution through QR decomposition)

#include "AK_Basic.h"

namespace MCMC{

/***** ****************************************************************************************** *****/
//
//  mean[dim]:         INPUT: whatsever
//                    OUTPUT: calculated mean of the normal approximation
//
//  tR[LT(dim)]:       INPUT: whatsever
//                    OUTPUT: lower triangular matrix t(R) stored in a packed format, where
//                            R is the part of the QR decomposition
//
//  log_det_R[1]:      INPUT: whatsever
//                    OUTPUT: log(det(R)) = log(prod R[j,j]) = sum log(abs(R[j,j]))
//                            * Dfull^{-1} = inverted variance of the normal approximation = t(R) %*% R
//                            * diagonal elements of R are not necessarily positive
//                              (number of negative diagonal elements is even)
//                              -> |Dfull|^{-1/2} = |R| = prod R[j, j]
//                              -> (-1/2) log|Dfull| = sum log(abs(R[j, j]))
//                         
//  QR[n + dim, dim]:  INPUT: whatsever
//                    OUTPUT: QR decomposition (is a packed format) of the working "design matrix"
//                            (in COLUMN major order)
//                            * R matrix is stored in the upper triangle of QR
//
//  uwork[n + dim]:    INPUT: whatsever
//                    OUTPUT: working "observations" for the LS solution
//
//  rsd[n + dim]:      INPUT: whatsever
//                    OUTPUT: residuals from the LS solution
//
//  tQu[n + dim]:      INPUT: whatsever
//                    OUTPUT: t(Q) %*% u, where u is a vector of working "observations"
//
//  rank[1]:           INPUT: whatsever
//                    OUTPUT: rank of the working "design" matrix,
//                            - value different from dim indicates problems
//
//  jpvt[dim]:         INPUT: whatsever
//                    OUTPUT: values possibly re-shuffled by 'dqrls'
//
//  QRaux[dim]:        INPUT: whatsever
//                    OUTPUT: auxiliary information needed to reconstruct the QR decomposition from the packed format
//
//  dwork[2 * dim]:    working array for 'dqrls'  
//
//  err[1]:            error flag
//
//  mean0[dim]:        initial guess for the mean of the normal approximation
//
//  uwork1[n]:         the first part of the working "observations" used in the LS solution
//                     (= standardized residuals from the GLMM model)
//
//  Zwork1[n, dim]:    the first part of the working "design matrix" used in the LS solution
//                     stored in COLUMN major order
//                     (= original GLMM design matrix modified by sqrt(var(y)) / phi)                   
// 
//  mu_prior[dim]:     prior mean
//
//  Li_prior[LT(dim)]: lower triangular matrix with decompositions of the prior inverse variance 
//                     (stored in a packed format in COLUMN major order)
//                     * D^{-1} = inverted prior variance = Li_prior %*% t(Li_prior)
//                     
//
//  n[1]:              number of observations
//
//  dim[1]:            dimension
//
//  half_factor[1]:    factor for step halving when using this routine for finding a mode in the context of GLMM
//                     * should be set to 1 for standard use 
//
//  caller:            name of the routine which has called this function 
//                     (used in error messages)
// 
/***** ****************************************************************************************** *****/

/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox_QR                                                             *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes: mean, decomposition (tR) of the inverted variance of the normal approximation
//  and log(det(tR)).
//
void
Moments_NormalApprox_QR(double* mean,
                        double* tR,
                        double* log_det_R,
                        double* QR,
                        double* uwork,
                        double* rsd,
                        double* tQu,
                        int*    rank,
                        int*    jpvt,
                        double* QRaux,
                        double* dwork,
                        int*    err,
                        const double* mean0,
                        const double* uwork1,
                        const double* Zwork1,
                        const double* mu_prior,
                        const double* Li_prior,
                        const int* n,
                        const int* dim,
                        const double* half_factor,
                        const char* caller);


/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox_QR (PROTOTYPE 2)                                               *****/
/***** ***************************************************************************************** *****/
//
//  This prototypes computes: log(det(tR)), where tR %*% R = inverted variance of the normal approximation
//
void
Moments_NormalApprox_QR(double* log_det_R,
                        double* QR,
                        int*    rank,
                        int*    jpvt,
                        double* QRaux,
                        double* dwork,
                        int*    err,
                        const double* Zwork1,
                        const double* Li_prior,
                        const int* n,
                        const int* dim,
                        const char* caller);


}    // end of namespace MCMC

#endif

