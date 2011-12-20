//
//  PURPOSE:   (Multivariate) GLMM, compute values of linear predictors and related quantities
//             from scratch
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   12/07/2009
//
//  FUNCTIONS:  
//     *   12/07/2009:  GLMM::linear_predictors  
//
//     *   27/02/2010:  GLMM::linear_predictors_fixed_updated
//         12/04/2010:    additional arguments added
//
//     *   07/08/2009:  GLMM::linear_predictor_fixed
//
//     *   07/08/2009:  GLMM::linear_predictor_random
//
//     *   07/08/2009:  GLMM::linear_predictor_zs
//
//     *   13/04/2010:  GLMM::linear_predictor_b_random_meanY
//
//     *   05/12/2011:  GLMM::linear_predictors_random_updated
//
// ======================================================================
//
#ifndef _GLMM_LINEAR_PREDICTORS_H_
#define _GLMM_LINEAR_PREDICTORS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"

#include "GLMM.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors                                                                   *****/
/***** ***************************************************************************************** *****/
//
//  eta_fixed[sum(n)]:      INPUT:  whatsever
//                         OUTPUT:  values of x[s,i,j]'beta (+ intercept)
//                                  filled by zeros if there are no betas/intercepts                               
//  
//  eta_random[sum(n)]:     INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'b[s,i] (+ random intercept)
//                                  filled by zeros if there are no random effects
//
//  eta[sum(n)]:            INPUT:  whatsever
//                         OUTPUT:  values of eta_fixed + eta_random
//
//  eta_zs[sum(n)]:         INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'shift_b[s]
//
//  N_s[R]:                 INPUT:  whatsever
//                         OUTPUT:  N_s[s] = number of observations for response s
//                                  sum(N_s) = sum(n)
//
//  N_i[I]:                 INPUT:  whatsever
//                         OUTPUT:  N_i[i] = number of observations for cluster i
//                                  sum(N_i) = sum(n)
//
//  OTHER VARIABLES:  see GLMM_MCMC.{h,cpp}
//
void
linear_predictors(double* eta_fixed,  
                  double* eta_random,      
                  double* eta,      
                  double* eta_zs,         
                  int*    N_s,
                  int*    N_i,
                  const double* X,    
                  const double* beta,      
                  const double* Z,  
                  const double* b,        
                  const double* shift_b,
                  const int*    p,       
                  const int*    fixedIntcpt,  
                  const int*    q,     
                  const int*    randIntcpt,  
                  const int*    n,       
                  const int*    R,            
                  const int*    I,     
                  const int*    dim_b,       
                  const int*    cumq_ri);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors_fixed_updated                                                     *****/
/***** ***************************************************************************************** *****/
//
//  eta_fixed[sum(n)]:      INPUT:  whatsever
//                         OUTPUT:  values of x[s,i,j]'beta (+ intercept)
//                                  filled by zeros if there are no betas/intercepts                               
//  
//  eta[sum(n)]:            INPUT:  whatsever
//                         OUTPUT:  values of eta_fixed + eta_random
//
//  meanY[sum(n)]:          INPUT:  whatsever
//                         OUTPUT:  updated values of E(Y | eta)
//
//  eta_random[sum(n)]:           linear predictor values derived from random effects
//
//  X[]:                          covariates for fixed effects
//                                see GLMM_MCMC.{h,cpp}
//
//  beta[sum(p + fixedIntcpt)]:   fixed effects
//
//  p[R]:                         number of covariates (intercept excluding) for each response
//
//  fixedIntcpt[R]:               0/1 indicating whether there is a fixed intercept for a specific response
//
//  dist[R]:                      types of response
//
//  n[I, R]:                      number of observations for each cluster and each response
//
//  R[1]:                         number of responses
//
//  I[1]:                         number of clusters
//
void
linear_predictors_fixed_updated(double* eta_fixed,  
                                double* eta,      
                                double* meanY,
                                const double* eta_random,
                                const double* X,    
                                const double* beta,      
                                const int*    p,       
                                const int*    fixedIntcpt,  
                                const int*    dist,
                                const int*    n,       
                                const int*    R,            
                                const int*    I);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictors_random_updated                                                    *****/
/***** ***************************************************************************************** *****/
//
//  eta_random[sum(n)]:     INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'b[s,i] (+ b[s,i,intercept])
//                                  filled by zeros if there are no b's/intercepts                               
//  
//  eta[sum(n)]:            INPUT:  whatsever
//                         OUTPUT:  values of eta_fixed + eta_random
//
//  meanY[sum(n)]:          INPUT:  whatsever
//                         OUTPUT:  updated values of E(Y | eta)
//
//  eta_fixed[sum(n)]:            linear predictor values derived from fixed effects
//
//  Z[]:                          covariates for random effects
//                                see GLMM_MCMC.{h,cpp}
//
//  b[I, sum(q + randIntcpt)]:    random effects
//
//  q[R]:                         number of random effects covariates (intercept excluding) for each response
//
//  randIntcpt[R]:                0/1 indicating whether there is a random intercept for a specific response
//
//  dist[R]:                      types of response
//
//  n[I, R]:                      number of observations for each cluster and each response
//
//  R[1]:                         number of responses
//
//  I[1]:                         number of clusters
//
//  dim_b[1]:                     dimension of random effects  = sum(q + randIntcpt)
//
void
linear_predictors_random_updated(double* eta_random,      
                                 double* eta,      
                                 double* meanY,
                                 const double* eta_fixed,  
                                 const double* Z,  
                                 const double* b,        
                                 const int*    q,     
                                 const int*    randIntcpt,  
                                 const int*    dist,
                                 const int*    n,       
                                 const int*    R,        
                                 const int*    I,
                                 const int*    dim_b);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_fixed                                                              *****/
/***** ***************************************************************************************** *****/
//
//  eta_fixed[sum(n), I]:   INPUT:  whatsever
//                         OUTPUT:  values of x[s,i,j]'beta (+ intercept)
//                                  filled by zeros if there are no betas/intercepts
//  
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictor_fixed function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_fixed(double* eta_fixed,
                       const double* X,    
                       const double* beta,
                       const int*    p,       
                       const int*    fixedIntcpt,
                       const int*    n,       
                       const int*    R,            
                       const int*    I);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_random                                                             *****/
/***** ***************************************************************************************** *****/
//
//  eta_random[sum(n), I]:  INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'b[s,i] (+ random intercept)
//                                  filled by zeros if there are no random effects
//
// 
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictor_random function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_random(double* eta_random,
                        const double* Z,     
                        const double* b,
                        const int*    q,        
                        const int*    randIntcpt,  
                        const int*    n,        
                        const int*    R,            
                        const int*    I,     
                        const int*    dim_b,       
                        const int*    cumq_ri);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_zs                                                                 *****/
/***** ***************************************************************************************** *****/
//  eta_zs[sum(n), I]:      INPUT:  whatsever
//                         OUTPUT:  values of z[s,i,j]'shift_b[s]
//
//  n[I]:                  number of observations for each longitudinal profile
//                         REMARK:  This is different from 'n' argument in linear_predictors function!!!
//                                  In linear_predictors_zs function, it is assumed that the number of observations
//                                  within each longitudinal profile is the same for all responses
//
//
//  OTHER VARIABLES:  see GLMM_longitClust.{h,cpp}
//
void
linear_predictor_zs(double* eta_zs, 
                    const double* Z,  
                    const double* shift_b,
                    const int*    q,     
                    const int*    randIntcpt,  
                    const int*    n,     
                    const int*    R,           
                    const int*    I,        
                    const int*    dim_b,       
                    const int*    cumq_ri);


/***** ***************************************************************************************** *****/
/***** GLMM::linear_predictor_b_random_meanY                                                     *****/
/***** ***************************************************************************************** *****/
//
//  It updates b, eta_random, eta, meanY under the assumption of the Gaussian distribution with identity link,
//  i.e., meanY = eta = eta_fixed + eta_random
//
//  It is primarily used in GLMM_updateRanEf_QR.
//
//  NOTATION:   dim_b = sum(q + randIntcpt)
//
//  b[dim_b]:              INPUT:  whatsever
//                        OUTPUT:  shift_b + scale_b * bscaled
//
//  eta_randomresp[R_c]:   INPUT:  pointers to eta_random of each response
//                                 eta_randomresp[s] is of length *nresp[s]
//                        OUTPUT:  each eta_randomresp[s] is updated using a new value of b
//                                 and shifted to point to the next cluster
//
//  etaresp[R_c]:          INPUT:  pointers to etaof each response
//                                 etaresp[s] is of length *nresp[s]
//                        OUTPUT:  each etaresp[s] is updated using a new value of b
//                                 and shifted to point to the next cluster
//
//  meanYresp[R_c]:        INPUT:  pointers to meanY of each response
//                                 meanYresp[s] is of length *nresp[s]
//                        OUTPUT:  each meanYresp[s] is updated using a new value of b
//                                 and is equal to eta_fixedresp[s] + eta_randomresp[s].
//                                 The pointer is also shifted to point to the next cluster
//
//  eta_fixedresp[R_c]:    pointers to eta_fixed for each response
//                         eta_fixedresp[s] is of length *nresp[s]
//                         ON EXIT: eta_fixedresp[s] is shifted to point to the next cluster
//
//  Zresp[R_c]:            pointers to Z matrix for each response
//                         ON EXIT: Zresp[s] is shifted to point to the next cluster
//
//  nresp[R_c]:            *nresp[s] determines the number of observations for the s-th response
//                         UNCHANGED ON EXIT
//
//  bscaled[dim_b]:        scaled values of random effects
//
//  shift_b[dim_b]:        shift vector for random effects
//
//  scale_b[dim_b]:        scale vector for random effects
//
//  q[R_c]:                
//
//  randIntcpt[R_c]:
//
//  R_c[1]:
// 
void
linear_predictor_gauss_b_random_meanY(double*  b,
                                      double** eta_randomresp,
                                      double** etaresp,
                                      double** meanYresp,
                                      double** eta_fixedresp,      // this is in fact const
                                      double** Zresp,              // this is in fact const
                                      int**    nresp,              // this is in fact const
                                      const double* bscaled, 
                                      const double* shift_b,
                                      const double* scale_b,
                                      const int*    q,
                                      const int*    randIntcpt,
                                      const int*    R_c);

}

#endif
