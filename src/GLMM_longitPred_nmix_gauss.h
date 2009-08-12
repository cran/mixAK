//
//  PURPOSE:   (Multivariate) GLMM, longitudinal prediction in the case
//             of all response variables being gaussian and
//             random effects having a normal mixture as distribution
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   07/08/2009
//
//  FUNCTIONS:  
//     * longitPred_nmix_gauss  07/08/2009:  start programming
//                              12/08/2009:  working version
//
// ======================================================================
//
#ifndef _GLMM_LONGITUDINAL_PREDICTION_NMIX_GAUSSIAN_H_
#define _GLMM_LONGITUDINAL_PREDICTION_NMIX_GAUSSIAN_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"
#include "AK_LAPACK.h"

#include "Dist_MVN.h"

#include "GLMM_linear_predictors.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::longitPred_nmix_gauss                                                               *****/
/***** ***************************************************************************************** *****/
//
// f_marg[sum(n)]:                 INPUT:   whatsever
//                                 OUTPUT:  f_marg[j, i] += f(y[0:(j-1), i])
//
// f_cond[sum(n)]:                 INPUT:   whatsever
//                                 OUTPUT:  f_cond[j, i] += f(y[0:(j-1), i] | hat(b)[0:(j-1), i])
// 
// f_ranef[sum(n)]:                INPUT:   whatsever
//                                 OUTPUT:  f_ranef[j, i] += f(hat(bscaled)[0:(j-1), i])
//                                          CORRECTION FOR SCALING (division by prod(scale_b)) NOT PERFORMED HERE
//
// eta_fixedresp[R_c]:             pointers to eta_fixed where eta_fixed for each response start
//
// eta_random[max(n)*R_c]:         working space to store eta_random for one longitudinal profile
//
// Y_crespP[R_c]:                  working pointers
//
// Y_drespP[1]:                    working pointers (not used in this function,
//                                 included only for prototype compatibility)
//
// eta_zsrespP[R_c]:               working pointers
//
// ZrespP[R_c]:                    working pointers
//
// err[1]:                         unchanged on exit if everything ok
//                                 non-zero if problems
//
// --------------- const arguments (some of them are not const to be able to use **) ------------------------
//
// Y_cresp[R_c]:                   pointers to Y_c where Y_c for each response starts
//
// Y_dresp[1]:                     NULL
//
// eta_zsresp[R_c]:                pointers to eta_zs where eta_zs for each response starts
//
// X[p_fi, sum(n), R_c]:           X matrices in ROW major order (intercept columns not included)
//   
// Zresp[R_c]:                     pointers to Z where Z matrix for each response starts
//
// SZitZiS[]:                      lower triangles of matrices S_s %*% t(Z_s[i,j]) %*% Z_s[i,j] %*% S_s, i=0,...,I-1, j=0,...,n[i]-1, s=0,...,R_c-1
//                                 where Z_s[i,j] are random effect design matrices (intercept columns included) for each prediction
//                                 and each response
//                                 !!! They are different from SZitZiS matrices used in GLMM::updateRanEf_nmix_gauss !!!
//
// ZiS[]:                          matrices Z_s[i,j] %*% S_s,
//                                 where Z_s[i,j] is the Z matrix for response s in the (i, j)-th prediction
//                                 (including possible column of ones for random intercept)
//                                 and S_s is the diagonal matrix with scale_b for the s-th response
//
// shift_b[dim_b]:                 shift for random effects
//
// scale_b[dim_b]:                 scale for random effects
//
// p[R_c]:                         number of columns in X matrix for each response
//
// fixedIntcpt[R_c]:               0/1 indicating whether a fixed intercept is included in the model for each response
//
// p_fi[R_c]:                      p_fi[i] = p[i] + fixedIntcpt[i]                        
//
// q[R_c]:                         number of columns in Z matrix for each response
//
// randIntcpt[R_c]:                0/1 indicating whether a random intercept is included in the model for each response
//
// q_ri[R_c]:                      q_ri[i] = q[i] + randIntcpt[i]
//
// cumq_ri[R_c]:                   cumsum(q_ri[i]) = cumsum(q[i] + randIntcpt[i])
//
// dim_b[1]:                       dimension of random effects
//
// LT_b[1]:                        length of lower triangle of random effect covariance matrices
//
// R_c[1]:                         number of continuous response variables
// 
// R_d[1]:                         NULL (number of discrete response variables -> it is assumed to be zero here)
//
// I[1]:                           number of longitudinal profiles
//
// n[I]:                           number of observations in each longitudinal profile
//                                 (it is assumed to be the same for each response type)
//
// max_n[1]:                       max(n)
//
// beta[p_fi]:                     value of fixed effects beta
//        
// sigma_eps[R_c]:                 value of residual standard deviations for each continuous response
// 
// K_b[1]:                         number of mixture components in the distribution of random effects
//
// w_b[K_b]:                       mixture weights for the distribution of random effects
//
// mu_b[dim_b, K_b]:               mixture means for the distribution of random effects
//
// Li_b[LT_b, K_b]:                Cholesky decompositions of mixture precision matrices for the distribution of random effects
//
void
longitPred_nmix_gauss(double* f_marg,            double* f_cond,            double* f_ranef,
                      double** eta_fixedresp,    double* eta_random,
                      double* log_dets_b,        double* dwork,             int* iwork,
                      double** Y_crespP,         int** Y_drespP,
                      double** eta_fixedrespP,   double** eta_zsrespP,
                      double** ZrespP,
                      int* err,
                      double** Y_cresp,          int** Y_dresp,
                      double** eta_zsresp,
                      const double* X,           double** Zresp,            const double* SZitZiS,    const double* ZiS,
                      const double* shift_b,     const double* scale_b,
                      const int* p,              const int* fixedIntcpt,
                      const int* q,              const int* randIntcpt,     const int* q_ri,          const int* cumq_ri,
                      const int* dim_b,          const int* LT_b,
                      const int* R_c,            const int* R_d,
                      const int* I,              const int* n,              const int* max_n,
                      const double* beta,        const double* sigma_eps,
                      const int* K_b,            const double* w_b,         const double* mu_b,        const double* Li_b);

}

#endif

