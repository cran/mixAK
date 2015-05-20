//
//  PURPOSE:   Main function to perform discriminant analysis of longitudinal profiles
//             based on models fitted with MCMC
//             - mixed responses allowed by this function
//             - also estimated values of random effects returned
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/04/2015
//
//  FUNCTIONS:  
//     * GLMM_longitDA2  13/04/2015:  Start working on it
//                       dd/mm/2015:  Working version available
//
// ======================================================================
//
//
#ifndef _GLMM_LONGITUDINAL_DISCRIMINANT_ANALYSIS_VERSION_TWO_H_
#define _GLMM_LONGITUDINAL_DISCRIMINANT_ANALYSIS_VERSION_TWO_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "NMix_Utils.h"
#include "NMix_updateAlloc.h"

#include "GLMM_create_ZS.h"
#include "GLMM_linear_predictors.h"
#include "GLMM_dY_meanY.h"
#include "GLMM_Deviance2.h"


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** GLMM_longitDA2                                                                            *****/
/***** ***************************************************************************************** *****/
//
//  Y_c
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
//  I[1]:
//
//  n[I, R]:                 number of observations in each longitudinal profile and each group 
//                           of repeated observations
//                           REMARK:  This is different to GLMM_longitDA, where it is assumed that 
//                                    within one longitudinal profile there is the same number of 
//                                    observations for each response variable!
//                                    Version used here is in agreement with all main fitting functions
//                                    (GLMM_MCMC, ...).
//                                    
//
//  X[]
//
//  p[R, nClust]:            numbers of columns of X matrices (intercept excluded) for each cluster and each response
// 
//  fixedIntcpt[R, nClust]:  0/1
//
//  q[R, nClust]:            numbers of columns of Z matrices (intercept excluded) for each cluster and each response
// 
//  randIntcpt[R, nClust]:  0/1
//
//  shiftScale_b[]
//
//  keepMCMC[nClust]:     length of MCMC for each cluster
//
//  keep_logf[1]          0 if logf_marg_i, logf_cond_i, logf_reff_i at each MCMC iteration are not kept, 
//                        != 0 if those are kept for each MCMC iteration
//
//  info[1]:              interval to print iteration info
//
//  Kmax_b[nClust]:       maximal number of K in the distribution of random effects for each cluster
//
//  f_marg[I, nClust]    estimated cluster-specific densities based on marginal approach
// 
//  f_cond[I, nClust]    estimated cluster-specific densities based on conditional approach
//
//  f_reff[I, nClust]    estimated cluster-specific densities based on random effect approach
//
//  bpred[I, dim_b_cl, nClust] estimated values of random effects (for each cluster - model)
//
//  zero_cond[I, nClust] number of MCMC iterations for which f_cond[i] in each cluster was equal to zero
//
//  zero_reff[I, nClust] number of MCMC iterayions for which f_reff[i] in each cluster was equal to zero
//

void
GLMM_longitDA2(const int*    nonSilent,
               double*       Y_c,                 /* it is in fact const, not const to be able to use ** */
               const int*    R_c,
               int*          Y_d,                 /* it is in fact const, not const to be able to use ** */
               const int*    R_d,
               const int*    dist,
               const int*    nClust,
               const int*    I,
               int*          n,                   /* it is in fact const, not const to be able to use ** */
               const double* X,
               const int*    p,
               const int*    fixedIntcpt,
               double*       Z,                   /* it is in fact const, not const to be able to use ** */
               const int*    q,
               const int*    randIntcpt,
               const double* shiftScale_b,
               const int*    distribution_b,
               const int*    keepMCMC,
	       const int*    keep_logf,
               const int*    info,
               const int*    Kmax_b,
               const double* chsigma_eps,
               const int*    chK_b,
               const double* chw_b,           
               const double* chmu_b,  
               const double* chLi_b,
               const double* chbeta,
               double*       f_marg,
               double*       f_cond,
               double*       f_reff,
               double*       logf_marg,
               double*       logf_cond,
               double*       logf_reff,
	       double*       bpred,
	       double*       logf_marg_i,
	       double*       logf_cond_i,
	       double*       logf_reff_i,
	       int*          nzero_marg,
	       int*          nzero_cond,
	       int*          nzero_reff,
               int*          err);

#ifdef __cplusplus
}
#endif

#endif
