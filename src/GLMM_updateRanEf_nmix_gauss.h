//
//  PURPOSE:   (Multivariate) GLMM, update of random effects in the case
//             of all response variables being gaussian and
//             random effects having a normal mixture as distribution
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//
//  FUNCTIONS:  
//     * updateRanEf_nmix_gauss  13/07/2009:  
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_RANDOM_EFFECTS_NMIX_GAUSSIAN_H_
#define _GLMM_UPDATE_RANDON_EFFECTS_NMIX_GAUSSIAN_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf_nmix_gauss                                                              *****/
/***** ***************************************************************************************** *****/
//
//  b[]:                   INPUT:  whatsever
//                        OUTPUT:  updated values of random effects
//
//  bscaled[]:             INPUT:  whatsever
//                        OUTPUT:  updated values of scaled random effects
//
//  eta_randomresp[R_c]:   INPUT:  pointers to places where eta_random for each response stars
//                        OUTPUT:  values to which it points are updated according to new values of random effects
//
//  eta_zsresp[R_c]:       INPUT:  pointers to places where eta_zs for each response stars
//                        OUTPUT:  unaltered
//
//  mu_full[dim_b]:        INPUT:  whatsever
//                        OUTPUT:  mean of the full conditional distribution of b[I-1]
//                                (also used as a working space during computation)
//
//  Li_full[LT(dim_b)]:    INPUT:  whatsever
//                        OUTPUT:  Cholesky decomposition of the precision matrix of the full conditional distribution of b[I-1]
//
//  log_dets[2]:           INPUT:  log_dets[0]: whatsever
//                                 log_dets[1]: -dim_b*log(sqrt(2pi))
//                        OUTPUT:  log_dets[0] = log(|Q_full[I-1]|^{1/2}) = sum(log(Li_full[I-1][j,j])),
//                                 where Q_full[I-1] is the precision matrix of the full conditional distribution of b[I-1]
//                                 log_dets[1]: unaltered
//
//  Qmu[dim_b*K]:          INPUT:  whatsever
//                        OUTPUT:  Q[k]*mu[k], k=0,...,K=1 for mixture components
//
//  dwork[dim_b]:          working array for Dist::rMVN2
//
//  Y_crespP[R_c]:         working array
//
//  Y_drespP:              NULL (it is here to get the same prototype for all updateRanEf functions)
//
//  eta_fixedrespP[R_c]:   working array
//
//  eta_randomrespP[R_c]:  working array
//
//  eta_zsrespP[R_c]:      working array
//
//  ZrespP[R_c]:           working array
//
//  nrespP[R_c]:           working array
//
//  err[1]                 INPUT:  whatsever
//                        OUTPUT:  unaltered if no problems, something > 0 if problems
//
//  Y_cresp[R_c]:         pointers to Y_c where Y_c for each response start
// 
//  Y_dresp:              NULL (it is here to get the same prototype for all updateRanEf functions)
//
//  eta_fixedresp[R_c]:   pointers to eta_fixedresp where eta_fixed for each response start      
//
//  Zresp[R_c]:           pointers to Z where Z matrix for each response start
//
//  SZitZiS[]:
//
//  shift[dim_b]:         shift for the random effects distribution
//
//  scale[dim_b]:         scale for the random effects distribution
//
//  q[R_c]:
//
//  randIntcpt[R_c]:
//
//  q_ri[R_c]:
//
//  cumq_ri[R_c]:
//
//  dim_b[1]:             dimension of random effects (= sum(q) + sum(randIntcpt))
//
//  LT_b[1]:              length of lower triangle of matrix dim_b x dim_b
//                        = (dim_b * (dim_b + 1)) / 2
//
//  R_c[1]
//
//  R_d[1]:               NULL
//
//  I[1]:
//
//  n[R_c]:               pointers to n where each response start
//
//  N_s[R_c]:
//
//  sigma[R_c]:           residual standard deviations for each response
//
//  K[1]:                 number of mixture components in the distribution of random effects
//  
//  mu[dim_b*K]:          mixture means for the distribution of random effects
//
//  Q[LT_b*K]:            mixture inverse variances for the distribution of random effects
//
//  r[I]:                 mixture component allocations for the distribution of random effects
//
void
updateRanEf_nmix_gauss(double* b,                 double* bscaled,           double** eta_randomresp,  
                       double* mu_full,           double* Li_full,           double* log_dets,     
                       double* Qmu,               double* dwork,
                       double** Y_crespP,         int** Y_drespP,        
                       double** eta_fixedrespP,   double** eta_randomrespP,  double** eta_zsrespP,
                       double** ZrespP,           int** nrespP,
                       int* err,
                       double** Y_cresp,          int** Y_dresp,
                       double** eta_fixedresp,    double** eta_zsresp, 
                       double** Zresp,            const double* SZitZiS,  
                       const double* shift,       const double* scale,
                       const int* q,              const int* randIntcpt,    const int* q_ri,      const int* cumq_ri,
                       const int* dim_b,          const int* LT_b,
                       const int* R_c,            const int* R_d,           const int* I,               
                       int** nresp,               const int* N_s,
                       const double* sigma,       
                       const int* K,              const double* mu,         const double* Q,      const int* r);

}  /** end of namespace GLMM **/

#endif
