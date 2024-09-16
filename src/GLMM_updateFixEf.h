//
//  PURPOSE:   (Multivariate) GLMM, update of fixed effects
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009 as GLMM_updateFixEf_gauss.h
//             21/10/2009 changed to GLMM_updateFixEf.h
//             26/10/2009 version which allows non-gaussian response as well
//                        (using Metropolis-Hastings with proposal obtained using one Newton-Raphson/Fisher scoring step)
//             29/03/2010 bug in shifting x_resp corrected
//
//  FUNCTIONS:  
//     * updateFixEf_gauss  11/07/2009:  CHANGED ON 21/10/2009 to updateFixEf
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_FIXED_EFFECTS_H_
#define _GLMM_UPDATE_FIXED_EFFECTS_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

#include "GLMM.h"
//#include "GLMM_fitted_Bernoulli_Logit.h"
//#include "GLMM_fitted_Poisson_Log.h"

#include "LogLik_Bernoulli_Logit.h"
#include "LogLik_Poisson_Log.h"

#include "MCMC.h"
#include "MCMC_Moments_NormalApprox.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateFixEf                                                                         *****/
/***** ***************************************************************************************** *****/
//
//  beta[sum(p_fi)]:           INPUT:  current values of beta
//                            OUTPUT:  updated values of beta
//
//  eta_fixed[sum(n)]:         INPUT:  current values of linear predictors based on fixed effects
//                            OUTPUT:  updated values of linear predictors based on fixed effects
//
//  eta[sum(n)]:               INPUT:  current values of linear predictors
//                            OUTPUT:  updated values of linear predictors
//
//  meanY[sum(n)]:             INPUT:  current values of response means
//                            OUTPUT:  updated response means
//
//  log_dets[2*(R_c + R_d)]:   INPUT:  log_dets[0,2,...]: whatsever
//                                     log_dets[1,3,...]: -p_fi[s]*log(sqrt(2pi))
//                            OUTPUT:  log_dets[0,2,...]: log(|Q[s]|^{1/2}) = sum(log(Li[s][j,j])),
//                                     where Q[s] is the precision matrix of the full conditional distribution of beta[s]
//                                     and Li[s] its Cholesky decomposition
//                                     log_dets[1,3,...]: unaltered
//
//  dwork[3*max(p_fi) + max(LT(p_fi)) + 2*max(N_s)]:  working array for 
//                             * Dist::rMVN2
//                             * (canonical) mean of the full conditional distribution
//                             * (Cholesky decomposition) of the precision matrix of the full cond. distr.
//                             * to store proposal beta (used only for discrete responses)
//                             * to store proposal eta_fixed (used only for discrete responses)
//                             * to store proposal mean_Y_d (used only for discrete responses)
//  
//  naccept[R_c + R_d]:        INPUT:  whatsever
//                            OUTPUT:  naccept[s] = naccept[s] from INPUT, if proposed value of beta[s] not accepted
//                                     naccept[s] = 1 + naccept[s] from INPUT, if proposed value of beta[s] accepted
//
//  err[1]                     INPUT:  whatsever
//                            OUTPUT:  unaltered if no problems, something > 0 if problems
//
//  Y_c[]:                    values of continuous response (sorted as argument Y_c of GLMM_MCMC function)
// 
//  Y_d[]:                    values of discrete response (sorted as argument Y_d of GLMM_MCMC function)
//
//  dY_d[]:                   additional double constants to evaluate the log-likelihood of discrete responses
//                            * this should have the same length as Y_d
//                            * for Bernoulli response:  currently not used
//                            * for Poisson response:    = log(y!)
//
//  eta_random[]:             values of random effect parts of linear predictors of continuous responses corresponding to Y
//
//  scale[]:                  vector (1, 1, ..., 1) having the same length as all fixed effects, i.e.,
//                            sum(p_fi)
//
//  X[]:                      design matrices for fixed effects (possible intercept not included)
//                            (sorted as argument X of GLMM_MCMC function)
//
//  XtX[]:                    FIRST PART (CORRESPONDING TO s = 0, ..., R_c - 1):
//                              lower triangles of matrices t(X_s) %*% X_s, where X_s is the design matrix 
//                              of the fixed effects (including possibly column of ones for a fixed intercept) 
//                              for response s = 0, ..., R_c - 1 (CONTINUOUS RESPONSES)
//
//                            SECOND PART (CORRESPONDING TO s = R_c, ..., R_c + R_d - 1):
//                              lower triangles of matrices x_s[i,j] %*% t(x_s[i,j]), i=0, ..., I-1, j = 0, ..., n[i] - 1,
//                              where x_s[i,j] is the (i,j)-th row (taken as column vector) of matrix X_s
//                              for responses s = R_c, ..., R_c + R_d - 1 (DISCRETE RESPONSES)
//
//  p[R_c + R_d]:             length of beta, intercept excluded, for each response
//
//  fixedIntcpt[R_c + R_d]:   0/1 indicating whether fixed intercept is included for particular response
//
//  p_fi[R_c + R_d]:          length of beta for each response = p + fixedInctp
//
//  R_c[1]:                   number of continuous responses
//
//  R_d[1]:                   number of discrete responses
//
//  dist[R_c + R_d]:          type of the distribution/link (see enum _GLMM_dist in GLMM.h)
//                            dist[0,...,R_c-1] is currently ignored as it is assumed that
//                            all continuous responses are gaussian with identity link
//
//  I[1]:                     number of clusters
// 
//  n[(R_c + R_d)*I]:         numbers of observations for response and each cluster (sorted as argument n of GLMM_MCMC function)
//
//  N_s[R_c + R_d]:           total number of observations for each response (sum(n) = sum(N_s))
//
//  sigma[R_c]:               residual standard deviations for each continuous response
//
//  Mbeta[sum(p_fi)]:         prior means of beta's
//
//  Pbeta[sum(p_fi)]:         prior precisions (inverse variances) of beta's
//
//  Pbeta_Mbeta[sum(p_fi)]:   prior precision TIMES prior mean for each beta
//
//  sqrt_tune_scale[R_d]:     square roots of the scale factor by which we multiply the proposal covariance matrix
//                            for each discrete response profile
//
//  log_sqrt_tune_scale[R_d]: log(sqrt_tune_scale)
//
void
updateFixEf(double* beta,              
            double* eta_fixed,       
            double* eta,  
            double* meanY,
            double* log_dets,            
            double* dwork,
            int*    naccept,
            int*    err,
            const double* Y_c,         
            const int*    Y_d,          
            const double* dY,
            const double* eta_random,  
            const double* scale,
            const double* X,         
            const double* XtX,
            const int*    p,              
            const int*    fixedIntcpt,  
            const int*    p_fi,      
            const int*    R_c,            
            const int*    R_d,          
            const int*    dist,
            const int*    I,              
            const int*    n,            
            const int*    N_s,
            const double* sigma,       
            const double* Mbeta,       
            const double* Pbeta,     
            const double* Pbeta_Mbeta,
            const double* sqrt_tune_scale,
            const double* log_sqrt_tune_scale);

}  /** end of namespace GLMM **/

#endif
