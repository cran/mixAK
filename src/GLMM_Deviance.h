//
//  PURPOSE:   Compute the GLMM marginal (random effects integrated out) log-likelihood
//             and conditional (given random effects) log-likelihood.
//             Marginal likelihood is calculated using the Laplacian approximation.
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   14/04/2010
//
//  FUNCTIONS:  
//     * GLMM_Deviance   ??/04/2010
//                       14/06/2010:  log_AK changed to log0_AK and exp_AK changed to exp0_AK
//                                    when calculating the log-likelihood
//                                    -> it leads to finite likelihood even with likelihood of 1e-305
//                       02/12/2010:  arguments sum_Yd_i and sum_Yd (typically containing sum(log(y!)) for Poisson response)
//                                    added to make corrections to log-likelihoods before exp to avoid exp(-Inf)
//                       02/12/2011:  argument iterate_to_mode added which allows for more than one Newton-Raphson step
//                                    when searching (for each k) for the mode of p(y_i | b, psi)*p(b | theta, u_i=k)
//
//     * TO DO:  
//       Commands:  *marg_ll_iP += *w_k * AK_Basic::exp0_AK(loglik_k);  (around line 239)
//                  *marg_ll_iP = AK_Basic::log0_AK(*marg_ll_iP); (around line 252)
//       In the case of the poisson response, it would be desirable to add log(y!) to loglik_k before exp-it    
//       and subtract it again after log the mixture sum. With higher counts, loglik_k is quite considerably negative
//       which leads to exp(-Inf) = 0 and then upon log again, we get -Inf log-likelihood
//     
//
// ======================================================================
//
#ifndef _GLMM_DEVIANCE_H_
#define _GLMM_DEVIANCE_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "AK_BLAS.h"

#include "MCMC_loglik_Zwork1_stres.h"
#include "MCMC_Moments_NormalApprox_QR.h"

namespace GLMM{

const int use_Hessian_in_bhat = 1;     // 0 --> in the Laplacian approximation, we use Hessian calculated
                                       //       in current value of b and not in the located mode bhat
                                       //       -> it is faster (Hessian is calculated only once, to locate bhat)
                                       // 1 --> in the Laplacian approximation, we use Hessian calculated
                                       //       in the located mode bhat
                                       //       -> Hessian must be calculated twice (once to locate bhat
                                       //          and once to calculate it in bhat)
// IMPORTANT NOTE:  When 'iterate_to_mode' is 1, value of use_Hessian_in_bhat does not have any effect
//                  and calculation proceeds as if use_Hessian_in_bhat = 1.

const int max_NRstep_Deviance = 10;            // see explanation of iterate_to_mode argument below
const double toler_NRstep_Deviance = 1e-5;     // see explanation of iterate_to_mode argument below
const int max_stephalf_Deviance = 10;         


/***** ********************************************************************** *****/
/***** GLMM:Deviance                                                          *****/
/***** ********************************************************************** *****/
//
// The function returns marginal (random effects integrated out) and conditional (given random effects)
// log-likelihood (overall and per grouped observations).
// The values on exit are NOT multiplied by -2 (to get the deviance).
//
/***** ********************************************************************** *****/
//
//  marg_ll[1]:        INPUT:  whatsever
//                    OUTPUT:  calculated value of the marginal log-likelihood
//                             = sum(marg_ll_i)
//
//  marg_ll_i[I]:      INPUT:  whatsever
//                    OUTPUT:  calculated value of the marginal log-likelihood for each cluster of grouped observations
//                             marg_ll_i[i] = log(sum(w[k] * marg_L_ik[i, k]))
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
//  cond_ll_i[I]:     INPUT:  whatsever
//                   OUTPUT:  calculated value of the conditional log-likelihood for each cluster of grouped observations
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
//  iterate_to_mode[1]:    0 --> Only one Newton-Raphson iteration is performed for each k when searching for the mode 
//                               of f(b) = p(y_i | b, psi)*p(b | theta, u_i=k).
//                        >0 --> At most GLMM::max_NRstep_Deviance Newton-Raphson steps are performed when searching for the mode of f(b).
//                               Iterations stops when |log f(b[t+1]) - log f(b[t])|/|log f(b[t+1])| < GLMM::toler_NRstep_Deviance.
//                                
//
void
Deviance(double* marg_ll,
         double* marg_ll_i,
         double* pi_ik,
         double* cond_ll,
         double* cond_ll_i,
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
         const int*    K,              
         const double* w,
         const double* logw,
         const double* mu,         
         const double* Li,
         const double* log_dets,
         const double* bscaled,
         const int*    iterate_to_mode);

}    // end of namespace GLMM

#endif
