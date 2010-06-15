//
//  PURPOSE:   (Multivariate) GLMM, update of random effects in the case
//             of random effects having a normal mixture as distribution
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   11/07/2009
//             29/03/2010  (validated in R, bug in shifting scale_resp pointer corrected)
//
//  FUNCTIONS:  
//     * updateRanEf_nmix_gauss                             13/07/2009
//     * updateRanEf created from updateRanEf_nmix_gauss on 26/10/2009
//
// ======================================================================
//
#ifndef _GLMM_UPDATE_RANDOM_EFFECTS_H_
#define _GLMM_UPDATE_RANDON_EFFECTS_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"

#include "GLMM.h"

#include "LogLik_Gauss_Identity.h"
#include "LogLik_Bernoulli_Logit.h"
#include "LogLik_Poisson_Log.h"

#include "MCMC.h"
#include "MCMC_Moments_NormalApprox.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::updateRanEf                                                                         *****/
/***** ***************************************************************************************** *****/
//
//  b[]:                         INPUT:  current values of random effects b
//                              OUTPUT:  updated values of random effects
//
//  bscaled[]:                   INPUT:  current values of scaled random effects bscaled
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
//                              OUTPUT:  log_dets_full[0] = log(|Q_full[I-1]|^{1/2}) = sum(log(Li_full[I-1][j,j])),
//                                       where Q_full[I-1] is the precision matrix of the full conditional distribution of b[I-1]
//                                       log_dets_full[1]: unaltered
//
//  dwork[K*dim_b + 5*dim_b + 3*LT_b + 2*max(N_i) + dim_b*dim_b]:  working array for     
//                         * Q[k] %*% mu[k], k=0,...,K-1 
//                         * Dist::rMVN2
//                         * (canonical) mean of the full conditional distribution
//                         * (canonical) mean of the reversal proposal distribution
//                         * (Cholesky decomposition) of the precision matrix of the full conditional distribution
//                         * (Cholesky decomposition) of the precision matrix of the reversal proposal distribution
//                         * information matrix given response
//                         * proposed values of bscaled
//                         * proposed values of b
//                         * backup of the full conditional (inverse) variance
//                           (stored in full matrix)
//                         * proposed values of eta_random
//                         * proposed values of mean_Y_d (for this, in fact, less than N_i slots are usually needed)
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
//  eta_zsrespP[R_c]:            working array (needed only for continuous responses)
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
//  eta_zsresp[R_c]:          pointers to places where eta_zs for each response starts 
//                            (it is needed only for continuous responses)
//
//  Zresp[R_c + R_d]:         pointers to Z where Z matrix for each response starts
//
//  SZitZiS[]:
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
//  cumq_ri[R_c + R_d]:
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
//  sigma[R_c]:               residual standard deviations for each response
//
//  K[1]:                     number of mixture components in the distribution of random effects
//  
//  mu[dim_b*K]:              mixture means for the distribution of random effects
//  
//  Q[LT_b*K]:                mixture inverse variances for the distribution of random effects
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
updateRanEf(double*  b,                 
            double*  bscaled,         
            double** eta_randomresp,  
            double** etaresp,    
	    double** meanYresp,
            double*  log_dets_full,     
            double*  dwork,
            double** Y_crespP,         
            int**    Y_drespP,     
            double** dYrespP,   
            double** eta_fixedrespP,   
            double** eta_randomrespP,    
            double** eta_zsrespP,
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
            double** eta_zsresp,                   // this is in fact const
            double** Zresp,                        // this is in fact const
            const double* SZitZiS,  
            const double* shift,       
            const double* scale,
            const int*    q,              
            const int*    randIntcpt,     
            const int*    q_ri,      
            const int*    cumq_ri,
            const int*    dim_b,          
            const int*    LT_b,
            const int*    R_c,            
            const int*    R_d,   
            const int*    dist,        
            const int*    I,               
            int**         nresp,                  // this is in fact const
            const int*    N_i,
            const double* sigma,       
            const int*    K,              
            const double* mu,         
            const double* Q,  
            const double* Li,    
            const double* log_dets,
            const int*    r,
            const double* sqrt_tune_scale,
            const double* log_sqrt_tune_scale);

}  /** end of namespace GLMM **/

#endif
