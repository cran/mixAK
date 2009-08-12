//
//  PURPOSE:   (Multivariate) GLMM, update of fixed effects in the case
//             of all response variables being gaussian
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//
//  FUNCTIONS:  
//     * updateFixEf_gauss  11/07/2009:  
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_FIXED_EFFECTS_GAUSSIAN_H_
#define _GLMM_UPDATE_FIXED_EFFECTS_GAUSSIAN_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateFixEf_gauss                                                                   *****/
/***** ***************************************************************************************** *****/
//
//  beta[sum(p_fi)]:           INPUT:  whatsever
//                            OUTPUT:  updated values of beta
//
//  eta_fixed[sum(n)]:         INPUT:  whatsever
//                            OUTPUT:  updated values of linear predictors based on fixed effects
//
//  mu_full[max(p_fi)]:        INPUT:  whatsever
//                            OUTPUT:  mean of the full conditional distribution of beta[R-1]
//                                     (also used as a working space during computation)
//
//  Li_full[LT(max(p_fi))]:    INPUT:  whatsever
//                            OUTPUT:  Cholesky decomposition of the precision matrix of the full conditional distribution of beta[R-1]
//                                     (also used as a working space during computation)
//
//  log_dets[2*R_c]:           INPUT:  log_dets[0,2,...]: whatsever
//                                     log_dets[1,3,...]: -p_fi[s]*log(sqrt(2pi))
//                            OUTPUT:  log_dets[0,2,...]: log(|Q[s]|^{1/2}) = sum(log(Li[s][j,j])),
//                                     where Q[s] is the precision matrix of the full conditional distribution of beta[s]
//                                     and Li[s] its Cholesky decomposition
//                                     log_dets[1,3,...]: unaltered
//
//  dwork[max(p_fi)]:         working array for Dist::rMVN2
//  
//  err[1]                     INPUT:  whatsever
//                            OUTPUT:  unaltered if no problems, something > 0 if problems
//
//  Y_c[]:                    values of continuous response (sorted as argument Y_c of GLMM_MCMC function)
// 
//  Y_d[]:                    NULL (it is here to get the same prototype for all updateFixEf functions)
//
//  eta_random[]:             values of random effect parts of linear predictors of continuous responses corresponding to Y
//
//  X[]:                      design matrices for fixed effects (possible intercept not included)
//                            (sorted as argument X of GLMM_MCMC function)
//
//  XtX[]:                    lower triangles of matrices t(X_s) %*% X_s, where X_s is the design matrix of the fixed effects 
//                            (including possibly column of ones for a fixed intercept) for response s (s=0,...,R-1)
//                            (sorted as argument XtX of GLMM_MCMC function)
//
//  p[R_c]:                   length of beta, intercept excluded, for each response
//
//  fixedIntcpt[R_c]:         0/1 indicating whether fixed intercept is included for particular response
//
//  p_fi[R_c]:                length of beta for each response = p + fixedInctp
//
//  R_c[1]:                   number of continuous responses
//
//  R_d[1]:                   NULL (it is here to get the same prototype for all updateFixEf functions)
//                            (R = R_c + R_d)
//
//  I[1]:                     number of clusters
// 
//  n[R_c*I]:                 numbers of observations for response and each cluster (sorted as argument n of GLMM_MCMC function)
//
//  N_s[R_c]:                 total number of observations for each response (sum(n) = sum(N))
//
//  sigma[R_c]:               residual standard deviations for each response
//
//  Pbeta[sum(p_fi)]:         prior precisions (inverse variances) of beta's
//
//  Pbeta_Mbeta[sum(p_fi)]:   prior precision TIMES prior mean for each beta
//
void
updateFixEf_gauss(double* beta,              double* eta_fixed,          
                  double* mu_full,           double* Li_full,         double* log_dets,            double* dwork,
                  int* err,
                  const double* Y_c,         const int* Y_d,      
                  const double* eta_random,  const double* X,         const double* XtX,
                  const int* p,              const int* fixedIntcpt,  const int* p_fi,      
                  const int* R_c,            const int* R_d,          const int* I,               
                  const int* n,              const int* N_s,
                  const double* sigma,       const double* Pbeta,     const double* Pbeta_Mbeta);

}  /** end of namespace GLMM **/

#endif
