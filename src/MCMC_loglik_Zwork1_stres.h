//
//  PURPOSE:   Subcode used primarily in GLMM_updateRanEf_QR.cpp
//             and in GLMM_Deviance.cpp
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:   20100412  created
//
//  FUNCTIONS:  
//     *   12/04/2010:  MCMC:loglik_Zwork1_stres (PROTOTYPE 1)
//     *   12/04/2010:  MCMC:loglik_Zwork1_stres (PROTOTYPE 2)
//     *   14/04/2010:  MCMC:loglik_Zwork1 (PROTOTYPE 2)
//     *   14/04/2010:  MCMC:loglik (PROTOTYPE 2)
//     *   25/01/2012:  MCMC:Zwork1_stres2UI 
//
// =================================================================================
//
#ifndef _MCMC_LOGLIK_ZWORKONE_STRES_H_
#define _MCMC_LOGLIK_ZWORKONE_STRES_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "GLMM.h"

#include "LogLik_Gauss_Identity.h"
#include "LogLik_Bernoulli_Logit.h"
#include "LogLik_Poisson_Log.h"

namespace MCMC{

/*** Calculate the upper part of the Z matrix and "observational" vector entering the LS solver   ***/
/*** Calculate the current value of the (conditional given random effects) likelihood             ***/
/*** ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ***/
//
//  NOTATION:    R     = *R_c + *R_d
//               N_i   = *nresp[0] + ... + *nresp[R-1]
//               dim_b = sum(q_ri)
//
//  loglik[1]:               INPUT: whatsever
//                          OUTPUT: calculated value of the log-likelihood (conditional given random effects)
//
//  b[dim_b]:                INPUT: whatsever
//                          OUTPUT: shifted and scaled values of random effects 
//                                  (PROTOTYPE 2)
//
//  Zwork1[N_i, dim_b]:      INPUT: whatsever
//                          OUTPUT: calculated working design matrix for the LS solver
//
//  stres[N_i]:             INPUT: whatsever
//                          OUTPUT: calculated working observations for the LS solver
//
//  sqrt_w_phi[N_i]:         INPUT: whatsever
//                          OUTPUT: outputs from LogLik functions
//
//  eta_randomresp[R]:      (PROTOTYPE 1)
//
//  meanYresp[R]:           (PROTOTYPE 1)
//
//  eta_random[N_i]:         INPUT: whatsever
//                          OUTPUT: new values of eta_random
//                                  (PROTOTYPE 2)
//
//  meanY[N_i]:              INPUT: whatsever
//                          OUTPUT: new values of meanY
//                                  (PROTOTYPE 2)
//
//  eta_fixedresp[R]:
//
//  dYresp[R]:
//
//  Y_cresp[R_c]:
//
//  Y_dresp[R_d]:
//
//  nresp[R]:
//
//  Zresp[R]:               (PROTOTYPE 2)
//
//  ZS[]:
// 
//  sigma[R_c]:
//
//  shift_b[dim_b]:         (PROTOTYPE 2)
//
//  scale_b[dim_b]:         (PROTOTYPE 2)
//
//  q_ri[R]:
//
//  dist[R]:
//
//  R_c[1]:
//
//  R_d[1]:
//
/***** ***************************************************************************************** *****/


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1_stres (PROTOTYPE 1)                                                  *****/
/***** ***************************************************************************************** *****/
//
//  This prototype uses supplied eta_random, eta_fixed, meanY to calculate likelihood etc.
//  - it does not modify eta_random, eta_fixed, meanY
// ==================================================================================================
void
loglik_Zwork1_stres(double*  loglik,
                    double*  Zwork1,
                    double*  stres,
                    double*  sqrt_w_phi,
                    int*     err,
                    double** eta_randomresp,             // this is in fact const
                    double** meanYresp,                  // this is in fact const
                    double** eta_fixedresp,              // this is in fact const
                    double** dYresp,                     // this is in fact const
                    double** Y_cresp,                    // this is in fact const
                    int**    Y_dresp,                    // this is in fact const
                    int**    nresp,                      // this is in fact const
                    const double* ZS,                    
                    const double* sigma,
                    const int* q_ri,
                    const int* dist,
                    const int* R_c,
                    const int* R_d);


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1_stres (PROTOTYPE 2)                                                  *****/
/***** ***************************************************************************************** *****/
//
//  This prototype updates also: b, eta_random, meanY
// ==================================================================================================
void
loglik_Zwork1_stres(double*  loglik,
                    double*  b,
                    double*  Zwork1,
                    double*  stres,
                    double*  sqrt_w_phi,
                    double*  eta_random,
                    double*  meanY,
                    int*     err,
                    double** eta_fixedresp,              // this is in fact const
                    double** dYresp,                     // this is in fact const
                    double** Y_cresp,                    // this is in fact const
                    int**    Y_dresp,                    // this is in fact const
                    int**    nresp,                      // this is in fact const
                    double** Zresp,                      // this is in fact const
                    const double* bscaled,
                    const double* ZS,                    
                    const double* sigma,
                    const double* shift_b,
                    const double* scale_b,
                    const int* q,
                    const int* randIntcpt,                          
                    const int* q_ri,
                    const int* dist,
                    const int* R_c,
                    const int* R_d);


/***** ***************************************************************************************** *****/
/***** MCMC::loglik_Zwork1 (PROTOTYPE 2)                                                         *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes only log-likelihood and Zwork1, stres is not computed.
//  It also updates b
// ==================================================================================================
void
loglik_Zwork1(double*  loglik,
              double*  b,
              double*  Zwork1,
              double*  sqrt_w_phi,
              int*     err,
              double** eta_fixedresp,              // this is in fact const
              double** dYresp,                     // this is in fact const
              double** Y_cresp,                    // this is in fact const
              int**    Y_dresp,                    // this is in fact const
              int**    nresp,                      // this is in fact const
              double** Zresp,                      // this is in fact const
              const double* bscaled,
              const double* ZS,                    
              const double* sigma,
              const double* shift_b,
              const double* scale_b,
              const int* q,
              const int* randIntcpt,                          
              const int* q_ri,
              const int* dist,
              const int* R_c,
              const int* R_d);


/***** ***************************************************************************************** *****/
/***** MCMC::loglik (PROTOTYPE 2)                                                                *****/
/***** ***************************************************************************************** *****/
//
//  This prototype computes only log-likelihood,  Zwork1 and stres are not computed.
//  It also updates b
// ==================================================================================================
void
loglik(double*  loglik,
       double*  b,
       int*     err,
       double** eta_fixedresp,              // this is in fact const
       double** dYresp,                     // this is in fact const
       double** Y_cresp,                    // this is in fact const
       int**    Y_dresp,                    // this is in fact const
       int**    nresp,                      // this is in fact const
       double** Zresp,                      // this is in fact const
       const double* bscaled,
       const double* sigma,
       const double* shift_b,
       const double* scale_b,
       const int* q,
       const int* randIntcpt,                          
       const int* q_ri,
       const int* dist,
       const int* R_c,
       const int* R_d);


/***** *************************************************** *****/
/***** MCMC::Zwork1_stres2UI                               *****/
/***** *************************************************** *****/
//
//  U[sum(q_ri)]:     INPUT:  whatsever
//                   OUTPUT:  calculated score vector
//
//  I[LT(sum(q_ri))]  INPUT:  whatsever
//                   OUTPUT:  lower triangle of the information matrix
//
void
Zwork1_stres2UI(double*  U,
                double*  I,
                int*     err,
                int**    nresp,                      // this is in fact const
                const double* Zwork1,
                const double* stres,
                const double* sqrt_w_phi,
                const double* ZS,
                const int* N_i,
                const int* q_ri,
                const int* dim_b,
                const int* dist,
                const int* R_c,
                const int* R_d);


}  // end of namespace MCMC

#endif

