//
//  PURPOSE:   Compute the GLMM marginal (random effects integrated out) log-likelihood
//             and conditional (given random effects) log-likelihood.
//             Marginal likelihood is calculated using the Laplacian approximation.
//
//             This is an updated version of the GLMM_Deviance.[h,cpp] where we 
//             also return estimated mode of the random effects distribution
//             for each set of grouped observations (predictions of individual random effects).
//
//             Differences to GLMM_Deviance
//                * new argument bhat (predictions of random effects)
//                * argument iterate_to_mode removed (it is always iterated towards the mode)
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/04/2015 by taking a copy of GLMM_Deviance.h
//
//  FUNCTIONS:  
//     * GLMM_Deviance2   20150414
//                       
//     * TO DO(?) (remains from GLMM_Deviance.[h,cpp]):  
//       Commands:  *marg_ll_iP += *w_k * AK_Basic::exp0_AK(loglik_k);  (around line 239)
//                  *marg_ll_iP = AK_Basic::log0_AK(*marg_ll_iP); (around line 252)
//       In the case of the poisson response, it would be desirable to add log(y!) to loglik_k before exp-it    
//       and subtract it again after log the mixture sum. With higher counts, loglik_k is quite considerably negative
//       which leads to exp(-Inf) = 0 and then upon log again, we get -Inf log-likelihood
//     
//
// ======================================================================
//
#ifndef _GLMM_DEVIANCE_TWO_H_
#define _GLMM_DEVIANCE_TWO_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"

#include "MCMC_loglik_Zwork1_stres.h"
#include "MCMC_Moments_NormalApprox_QR.h"

#include "NMix.h"
#include "Dist_MVT.h"
#include "Dist_mixMVN.h"


namespace GLMM{

const int max_NRstep_Deviance2 = 10;            // see explanation of iterate_to_mode argument below
const double toler_NRstep_Deviance2 = 1e-5;     // see explanation of iterate_to_mode argument below
const int max_stephalf_Deviance2 = 10;         


/***** ********************************************************************** *****/
/***** GLMM:Deviance2                                                          *****/
/***** ********************************************************************** *****/
//
// The function returns marginal (random effects integrated out),
// conditional (given random effects) log-likelihood (overall and per grouped observations)
// and calculated predictions (approx. mode of the appropriate distribution)
// of random effects.
//
// The values on exit are NOT multiplied by -2 (to get the deviance).
//
/***** ********************************************************************** *****/
//
//  marg_ll[1]:        INPUT:  whatsever
//                    OUTPUT:  calculated value of the marginal log-likelihood
//                             = sum(marg_ll_i)
//
//  marg_L_i[I]:      INPUT:  whatsever
//                    OUTPUT:  calculated value of the marginal likelihood (not log!!!)
//                             for each set of grouped observations
//                             marg_L_i[i] = sum(w[k] * marg_L_ik[i, k])
//
//  pi_ik[K, I]:       INPUT:  whatsever
//                    OUTPUT:  calculated value of the marginal likelihood for each cluster of grouped observations
//                             given the mixture component times mixture weight, i.e.,
//                             pi_ik[k, i] = w[k] * marg_L_ik[K, i]
//
//  cond_ll[1]:       INPUT:  whatsever
//                   OUTPUT:  calculated value of the conditional log-likelihood
//                            = sum(cond_ll_i)
//
//  cond_L_i[I]:     INPUT:  whatsever
//                   OUTPUT:  calculated value of the conditional likelihood (not log!!)
//                            for each set of grouped observations
//                            given bpred_i (THIS IS DIFFERENT FROM WHAT RETURNS GLMM_Deviance which
//                            uses supplied bscaled)
//
//  bpred_i[I, dim_b]  INPUT:  whatsever
//                    OUTPUT:  calculated predictions of random effects
//
//  bpredscaled_i[I, dim_b]  INPUT:  STARTING values for algorithm that looks for the modes,
//                                   eta's etc. should correspond to it
//                          OUTPUT:  calculated predictions of random effects (scaled ones)
//
//  reff_ll[1]:       INPUT:  whatsever
//                   OUTPUT: calculated value of the random-effects log-likelihood evaluated in bhat
//
//  reff_L_i[I]:     INPUT:  whatsever
//                   OUTPUT:  calculated value of the random-effects likelihood (not log!!) evaluated in bpred_i
//                            for each set of grouped observations
//
//  nzero_marg_i[I]:  INPUT:  whatsever (typically number of zero L_marg's in previous MCMC iterations)
//                   OUTPUT: +1 on those places where marg_L_i[i] is zero in which case on exit marg_ll_i[i] = AK_Basic::_LOG_ZERO0
//
//  nzero_cond_i[I]:  INPUT:  whatsever (typically number of zero L_cond's in previous MCMC iterations)
//                   OUTPUT: +1 on those places where cond_L_i[i] is zero in which case on exit cond_ll_i[i] = AK_Basic::_LOG_ZERO0
//
//  nzero_reff_i[I]:  INPUT:  whatsever (typically number of zero L_reff's in previous MCMC iterations)
//                   OUTPUT: +1 on those places where reff_L_i[i] is zero in which case on exit reff_ll_i[i] = AK_Basic::_LOG_ZERO0
//
//  stres[sum(N_i)]:  INPUT:  whatsever
//                   OUTPUT:  standardized residuals for each observation
//                            cluster index is the major index to loop,
//                            the second index to loop is the response profile index
//
//  sqrt_w_phi[sum(N_i)]:  INPUT:  whatsever
//                        OUTPUT:  sqrt(var(Y_{i,s,j}|theta, b)) / phi_s, where phi_s is the dispersion parameter of the s-th response
//                                 cluster index (i) is the major index to loop,
//                                 the second index to loop is the response profile index (s)
//
//
//  iterNum[1]         Iteraion number (for error messages)
//
//
void
Deviance2(double* marg_ll,
	  double* marg_ll_i,
          double* marg_L_i,
          double* pi_ik,
          double* cond_ll,
          double* cond_ll_i,
          double* cond_L_i,
          double* bpred_i,
	  double* bpredscaled_i,
          double* reff_ll,
          double* reff_ll_i,
          double* reff_L_i,
	  int*    nzero_marg_i,
          int*    nzero_cond_i,
          int*    nzero_reff_i,
          double* stres,
          double* sqrt_w_phi,
          double** Y_crespP,
          int**    Y_drespP,
          double** dYrespP,   
          double** eta_fixedrespP,
          double** eta_randomrespP,
          double** meanYrespP,
          double** ZrespP,           
          int**    nrespP,
          int*     iwork,
          double*  dwork,
          int*     err,
          double** Y_cresp,              // this is in fact const
          int**    Y_dresp,              // this is in fact const
          double** dYresp,               // this is in fact const
          double** eta_fixedresp,        // this is in fact const
          double** eta_randomresp,       // this is in fact const
          double** meanYresp,            // this is in fact const
          double** Zresp,                // this is in fact const
          int**    nresp,                // this is in fact const
          const double* ZS,
          const double* shift,       
          const double* scale,
	  const double* scaleProd,
	  const double* logscaleSum,
          const int*    q,              
          const int*    randIntcpt,     
          const int*    q_ri,      
          const int*    dim_b,          
          const int*    LT_b,
          const int*    R_c,            
          const int*    R_d,   
          const int*    dist,        
          const int*    I,               
          const int*    N_i,
          const int*    max_N_i,
          const int*    l_ZS,
          const double* sigma,
          const int*    distribution_b,
          const int*    K,              
          const double* w,
          const double* logw,
          const double* mu,         
          const double* Li,
          const double* Q,                 // only needed for a mixture of MVT in the random effects distribution
          const double* df,                // only needed for a mixture of MVT in the random effects distribution
          const double* log_dets,
          const int*    iterNum);

}    // end of namespace GLMM

#endif
