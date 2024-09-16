//
//  PURPOSE:   Main function to perform discriminant analysis of longitudinal profiles
//             based on models fitted with MCMC
//             - only continuous responses allowed by this function
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/08/2009 as GLMM_longitClust.h
//             28/10/2009: renamed to GLMM_longitDA.h
//
//  FUNCTIONS:  
//     * GLMM_longitDA  06/08/2009:  Start working on it
//                      12/08/2009:  Version for all responses being continuous working (named as GLMM_longitClust)
//                      28/10/2009:  Renamed to GLMM_longitDA
//                      02/11/2009:  Matrices S %*% t(Zi) %*% Zi %*% S are computed directly in C++ code
//
// ======================================================================
//
//
#ifndef _GLMM_LONGITUDINAL_DISCRIMINANT_ANALYSIS_H_
#define _GLMM_LONGITUDINAL_DISCRIMINANT_ANALYSIS_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Error.h>
#include <R_ext/RS.h>

#include "AK_Basic.h"

#include "GLMM_longitPred_nmix_gauss.h"
#include "GLMM_linear_predictors.h"
#include "GLMM_create_SZitZiS_4longitDA.h"
#include "GLMM_scale_ZitZi.h"
#include "GLMM_create_ZiS.h"

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_longitDA                                                                             *****/
/***** ***************************************************************************************** *****/
//
//  Yc
//
//  R_c
//
//  Y_d
//
//  R_d
//
//  dist[R]
//
//  nClust[1]:               number of clusters
//
//  I[1]
//
//  n[I]:                    number of observations in each longitudinal profile
//                           REMARK:  In contrast to GLMM_MCMC, it is assumed that within one longitudinal profile
//                                    there is the same number of observations for each response variable!
//                                    
//
//  X[]
//
//  p[R, nClust]:            numbers of columns of X matrices (intercept excluded) for each cluster and each response
// 
//  fixedIntcpt[R, nClust]:  0/1
//
//  SZitZiS[]:               INPUT:  lower triangles of matrices t(Z_s[i]) %*% Z_s[i],
//                                   where Z_s[i] is the design matrix of the random effects
//                                   (including possibly column of ones for a random intercept), sorted in this way:
//                                   t(Z_0[0]) %*% Z_0[0], ..., t(Z_{R-1}[0]) %*% Z_{R-1}[0],
//                                   ...
//                                   t(Z_0[I-1]) %*% Z_0[I-1], ..., t(Z_{R-1}[I-1]) %*% Z_{R-1}[I-1], 
//                                   !!!! It is somewhat different from the same argument of GLMM_MCMC.
//                                   !!!! It contains such matrices for all instances of prediction.
//
//                           OUTPUT: lower triangles of matrices S %*% t(Z_s[i]) %*% Z_s[i] %*% S,
//                                   where S is the diagonal matrix with scale_b on a diagonal  
//
//         ARGUMENT SZitZiS HAS BEEN REMOVED ON 02/11/2009
//
//
//  q[R, nClust]:            numbers of columns of Z matrices (intercept excluded) for each cluster and each response
// 
//  randIntcpt[R, nClust]:  0/1
//
//  shiftScale_b[]
//
//  keepMCMC[nClust]:     length of MCMC for each cluster
//
//  info[1]:              interval to print iteration info
//
//  Kmax_b[nClust]:       maximal number of K in the distribution of random effects for each cluster
//
//  pi_marg[I, nClust]    estimated cluster probabilities based on marginal approach
// 
//  pi_cond[I, nClust]    estimated cluster probabilities based on conditional approach
//
//  pi_ranef[I, nClust]   estimated cluster probabilities based on random effect approach
//
void
GLMM_longitDA(double* Y_c,                       /* it is in fact const, not const to be able to use ** */
              const int* R_c,
              int* Y_d,                          /* it is in fact const, not const to be able to use ** */
              const int* R_d,
              const int* dist,
              const int* nClust,
              const int* I,
              const int* n,
              const double* X,
              const int* p,
              const int* fixedIntcpt,
              double* Z,                         /* it is in fact const, not const to be able to use ** */
              const int* q,
              const int* randIntcpt,
              const double* shiftScale_b,
              const int* keepMCMC,
              const int* info,
              const int* Kmax_b,
              const double* chsigma_eps,
              const int* chK_b,
              const double* chw_b,           
              const double* chmu_b,  
              const double* chLi_b,
              const double* chbeta,
              double* pi_marg,
              double* pi_cond,
              double* pi_ranef,
              int* err);

#ifdef __cplusplus
}
#endif

#endif
