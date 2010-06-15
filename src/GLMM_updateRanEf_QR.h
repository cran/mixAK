//
//  PURPOSE:   (Multivariate) GLMM, update of random effects in the case
//             of random effects having a normal mixture as distribution
//             * change compared to GLMM_updateRanEf:
//               moments of the proposal distribution are found through QR decomposition
//               and solution of the corresponding least-squares problem
// 
//             * it is also not necessary to make such a strict distinction between normal and non-normal response
//
//             * compute also (if asked for) the GLMM (marginal with random effects integrated out) 
//               log-likelihoods (for each grouped observations and full model) obtained
//               using the Laplacian approximation to the integral over random effects
//
//             * the marginal log-likelihood is defined as log L_i(beta),
//               where L_i(beta) = \int\prod_{r=1}^R\prod_{j=1}^{n_{i,r}} p(y_{i,r,j} | beta, b_i) p(b_i) db_i
//                               = \sum_{k=1}^K w_k \int\prod_{r=1}^R\prod_{j=1}^{n_{i,r}} p(y_{i,r,j} | beta, b_i) N(b_i | mu_k, D_k) db_i
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/04/2010
//
//  FUNCTIONS:  
//     * updateRanEf_QR   13/04/2010
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_RANDOM_EFFECTS_QR_H_
#define _GLMM_UPDATE_RANDON_EFFECTS_QR_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

#include "GLMM_linear_predictors.h"
#include "GLMM_copy_shift_eta_meanY_Zresp.h"

#include "MCMC.h"
#include "MCMC_loglik_Zwork1_stres.h"
#include "MCMC_Moments_NormalApprox_QR.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf_QR                                                                      *****/
/***** ***************************************************************************************** *****/
//
//  b[]:                         INPUT:  current values of random effects b
//                              OUTPUT:  updated values of random effects
//
//  bscaled[]:                   INPUT:  current values of scaled random effects
//                              OUTPUT:  updated values of scaled random effects
//
//  eta_randomresp[R_c + R_d]:   INPUT:  pointers to places where eta_random for each response starts
//                              OUTPUT:  values to which it points are updated according to new values of random effects
//
//  etaresp[R_c + R_d]:          INPUT:  pointers to places where eta for each response starts
//                              OUTPUT:  values to which it points are updated according to new values of random effects
//
//  meanYresp[R_c + R_d]:        INPUT:  pointers to places where meanY for each response starts
//                              OUTPUT:  values to which it points are updated according to new values of random effects
//
//  log_dets_full[2]:            INPUT:  log_dets_full[0]: whatsever
//                                       log_dets_full[1]: -dim_b*log(sqrt(2pi))
//                              OUTPUT:  log_dets_full[0] = log(|Q_full[I-1]|^{1/2}) = sum(log(abs(tR_full[I-1][j,j]))),
//                                       where Q_full[I-1] is the precision matrix of the full conditional distribution of b[I-1]
//                                       log_dets_full[1]: unaltered
//
//  iwork[dim_b]:                working array used as jpvt argument in MCMC::Moments_NormalApprox
//
//  dwork[]:                     working array
//
//  Y_crespP[R_c]:               working array
//
//  Y_drespP[R_d]:               working array           
//
//  dYrespP[R_c + R_d]:          working array           
//
//  eta_fixedrespP[R_c + R_d]:   working array
//
//  eta_randomrespP[R_c + R_d]:  working array
//
//  etarespP[R_c + R_d]:         working array
//
//  meanYrespP[R_c + R_d]:       working array           
//
//  ZrespP[R_c + R_d]:           working array
//
//  nrespP[R_c + R_d]:           working array
//
//  naccept[I]:            INPUT:  whatsever
//                        OUTPUT:  naccept[i] = naccept[i] from INPUT, if proposed value of b[i] not accepted
//                                 naccept[i] = 1 + naccept[i] from INPUT, if proposed value of b[i] accepted
//                              
//  err[1]                 INPUT:  whatsever
//                        OUTPUT:  unaltered if no problems, something > 0 if problems
//
//  Y_cresp[R_c]:             pointers to Y_c where Y_c for each response starts
// 
//  Y_dresp[R_d]:             pointers to Y_d where Y_d for each response starts
//
//  dYresp[R_c + R_d]:        pointers to dY where dY for each response starts
//
//  eta_fixedresp[R_c + R_d]: pointers to eta_fixedresp where eta_fixed for each response starts      
//
//  Zresp[R_c + R_d]:         pointers to Z where Z matrix for each response starts
//
//  ZS[]:                     Z[i,s] %*% S[s] (i < I, s< R_c + R_d) matrices
//
//  shift[dim_b]:             shift for the random effects distribution
//
//  scale[dim_b]:             scale for the random effects distribution
//
//  q[R_c + R_d]:
//
//  randIntcpt[R_c + R_d]:
//
//  q_ri[R_c + R_d]:
//
//  dim_b[1]:                 dimension of random effects (= sum(q) + sum(randIntcpt))
//
//  LT_b[1]:                  length of lower triangle of matrix dim_b x dim_b
//                            = (dim_b * (dim_b + 1)) / 2
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
//  nresp[R_c + R_d]:         pointers to n where each response start
//
//  N_i[I]:                   total number (for all response types) of observations per cluster
//
//  max_N_i[1]:               max(N_i)
//
//  l_ZS[I]:                  length of {Z[i,s] %*% S[s], s < R_c + R_d} for each i < I
//
//  sigma[R_c]:               residual standard deviations for each continuous response
//
//  //K[1]:                     number of mixture components in the distribution of random effects
//  
//  mu[dim_b*K]:              mixture means for the distribution of random effects
//  
//  Li[LT_b*K]:               Cholesky decompositions of mixture inverse variances for the distribution of random effects
//
//  log_dets[2*K]:            log_dets based on Q and Li matrices
//
//  r[I]:                     mixture component allocations for the distribution of random effects
//
//  sqrt_tune_scale[1]:       square root of the scale factor by which we multiply the proposal covariance matrix
//                            when there are some discrete response profile
//
//  log_sqrt_tune_scale[1]:   log(sqrt_tune_scale)
//
// ***************************************************************************************************************************
void
updateRanEf_QR(double* b,
               double* bscaled,
               double** eta_randomresp,  
               double** etaresp,
	       double** meanYresp,
               double*  log_dets_full,
               int*     iwork,     
               double*  dwork,
               double** Y_crespP,         
               int**    Y_drespP,     
               double** dYrespP,   
               double** eta_fixedrespP,   
               double** eta_randomrespP,
               double** etarespP,
               double** meanYrespP,    
               double** ZrespP,           
               int**    nrespP,
               int*     naccept,
               int*     err,
               double** Y_cresp,                      // this is in fact const
               int**    Y_dresp,                      // this is in fact const
               double** dYresp,                       // this is in fact const
               double** eta_fixedresp,                // this is in fact const
               double** Zresp,                        // this is in fact const
               const double* ZS,        
               const double* shift,       
               const double* scale,
               const int*    q,              
               const int*    randIntcpt,     
               const int*    q_ri,      
               const int*    dim_b,          
               const int*    LT_b,
               const int*    R_c,            
               const int*    R_d,   
               const int*    dist,        
               const int*    I,               
               int**         nresp,                  // this is in fact const           
               const int*    N_i,
               const int*    max_N_i,
               const int*    l_ZS,
               const double* sigma,       
               const double* mu,         
               const double* Li,
               const double* log_dets,
               const int*    r,
               const double* sqrt_tune_scale,
               const double* log_sqrt_tune_scale);

}    // end of namespace GLMM

#endif 
