//
//  PURPOSE:   Normal mixture model, smaller utilities
//
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:      20071126  created
//            20120126  stuff needed for MVT mixture added
//            20150327  modifications to allow for situation 
//                      when a factor covariate influences the mixture weights
//
//  FUNCTIONS:  
//     * w2logw                       26/11/2007:  
//
//     * Li2log_dets                  10/02/2010:
//
//     * wLi2w_dets                   08/11/2008:
//
//     * Li2Q                         26/11/2007:  
//
//     * Li2Sigma                     27/11/2007:  
//
//     * Li2sigma                     09/11/2008:
//
//     * muLi2beta_sigmaR2            09/11/2008:
//
//     * Moments                      27/11/2007:
//
//     * ySum_j                       26/11/2007:     
//                                    19/12/2007:  REMOVED
//     * ySumBar_j                    19/12/2007:  CREATED FROM ySum_j
//
//     * SS_j                         26/11/2007:  
//
//     * ySum_SSm_j                   13/02/2008:
//
//     * prior_derived                08/07/2009:
//
//     * init_derived                 08/07/2009:
//
// ====================================================================================================
//
#ifndef _NMIX_UTILS_H_
#define _NMIX_UTILS_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"
#include "Stat_BLA.h"
#include "Dist_Wishart.h"

#include "NMix.h"

namespace NMix{


/***** ***************************************************************************************** *****/
/***** NMix::w2logw                                                                              *****/
/***** ***************************************************************************************** *****/
void
w2logw(double* logw,  const double* w,  const int* K,  const int* nxw);


/***** ***************************************************************************************** *****/
/***** NMix::Li2log_dets                                                                         *****/
/***** ***************************************************************************************** *****/
//
//  log_dets[2, K]:   INPUT:   whatsever
//                    OUTPUT:  log_dets[0] = log(|Sigma[0]|^{-1/2}) = sum(log(Li_{0}[i,i]))
//                             log_dets[2] = log(|Sigma[1]|^{-1/2}) = sum(log(Li_{1}[i,i]))
//                             ...
//                             log_dets[1], log_dets[3], ... is unaltered
//
//  Li[LT(p), K]
//
//  K[1]
//
//  p[1]
//
void
Li2log_dets(double* log_dets,  const double* Li,  const int* K,  const int* p);


/***** ************************************************************************************ *****/
/***** NMix::wLi2w_dets                                                                     *****/
/***** ************************************************************************************ *****/
//
// w_dets[nxw * K]   INPUT:  whatsever
//                  OUTPUT: w_dets[i, k] = (2*pi)^(-p/2) * |Li[k]| * w[i, k]
//
// w[nxw * K]        mixture weights
//
// Li[K*LT(p)] Cholesky decompositions of mixture precision matrices
//
// K[1]        number of mixture components
//
// p[1]        dimension
//
// nxw[1]      number of factor covariates on weights
//
void
wLi2w_dets(double* w_dets,  const double* w,  const double* Li,  const int* K,  const int* p, const int* nxw);


/***** ***************************************************************************************** *****/
/***** NMix::Li2Q                                                                                *****/
/***** ***************************************************************************************** *****/
//
// Q[j] = Li[j] %*% t(Li[j])
//
// Q[LT(p), K]   OUTPUT:  Symmetric matrices stored in a packed form
//
// Li[LT(p), K]  INPUT:   Lower triangular matrices stored in a packed form
//
void
Li2Q(double* Q,  const double* Li,  const int* K, const int* p);


/***** ***************************************************************************************** *****/
/***** NMix::Li2Sigma                                                                            *****/
/***** ***************************************************************************************** *****/
//
//  Sigma[j] = (Li[j] %*% t(Li[j]))^{-1}
//
//  Sigma[LT(p), K]   OUTPUT:  Symmetric matrices stored in a packed form
//
//  Li[LT(p), K]      INPUT:   Lower triangular matrices stored in a packed form
//
void
Li2Sigma(double* Sigma,  int* err,  const double* Li,  const int* K,  const int* p);


/***** ************************************************************************************ *****/
/***** NMix::Li2sigma                                                                       *****/
/***** ***************** ****************************************************************** *****/
//
// sigma[j] = 1 / Li[j] 
//
// sigma[K]    OUTPUT:  standard deviations
//
// Li[K]       INPUT:   inverted standard deviations
//
void
Li2sigma(double* sigma,  const double* Li,  const int* K);


/***** ************************************************************************************ *****/
/***** NMix::muLi2beta_sigmaR2                                                              *****/
/***** ***************** ****************************************************************** *****/
//
// beta[K * p x p]   OUTPUT:  coefficients of full conditional distributions for all components
//
// sigmaR2[K * p]    OUTPUT:  residual variances (of full conditional distributions) for all components
//
// work[LT(p) + LT(p-1)]
//
void
muLi2beta_sigmaR2(double* beta,  double* sigmaR2,   double* work,
                  const int* K,  const double* mu,  const double* Li,  
                  const int* p,  const int* p_p,    const int* LTp);


/***** ***************************************************************************************** *****/
/***** NMix::Moments                                                                             *****/
/***** ***************************************************************************************** *****/
//
// Compute first two moments, standard deviations and correlations of a normal mixture
// and of a normal mixture for original data 
//
// * Mean = sum_j w_j * mu_j
// * Var  = sum_j w_j (Sigma_j + (mu_j - Mean)*t(mu_j - Mean)')
//
// For MVT distribution: 
//          Mean of one mixture component is equal to mu_j provided that df_j > 1.
//          This function calculates "Mean" as sum_j w_j * mu_j even if df_j <= 1.
//
//          Variance of one mixture component is equal to (df_j / (df_j - 2)) * Sigma_j provided that df_j > 2.
//          If df_j <= 2, this function calculates "Var" using 2.001/0.001 * Sigma_j.
//
//          Hence for MVT, use calculated results with care!
//
// Mean[p, nxw]
//
// Var[LT(p), nxw]
//
// Corr[LT(p), nxw]    Diagonal:      standard deviations = sqrt(Var[i,i])
//                  Off-diagonal:  correlations = Var[i,l]/sqrt(Var[i,i]*Var[l,l])
//
// MeanData[p, nxw]        MeanData = shift + scale*Mean
//
// VarData[LT(p), nxw]     VarData = diag(scale) %*% Var %*% diag(scale)
//
// CorrData[LT(p), nxw]    Diagonal:      standard deviations = sqrt(VarData[i,i])
//                      Off-diagonal:  correlations = VarData[i,l]/sqrt(VarData[i,i]*VarData[l,l]) = Var[i,l]/sqrt(Var[i,i]*Var[l,l])
//
// w[K, nxw]
//
// mu[p, K]
//
// Sigma[LT(p), K]
//
// K[1]
//
// shift[p]  
//
// scale[p]
//
// p[1]
//
// nxw[1]
//
void
Moments(double* Mean,         
        double* Var,          
        double* Corr,
        double* MeanData,     
        double* VarData,      
        double* CorrData,
        const int* distribution,
        const double* w,      
        const double* mu,     
        const double* Sigma,  
        const double* df,
        const int*    K,  
        const double* shift,  
        const double* scale,  
        const int* p,
	const int* nxw);


/***** ***************************************************************************************** *****/
/***** NMix::ySumBar_j                                                                           *****/
/***** ***************************************************************************************** *****/
//
// mixsumy[p, K]
//
// mixbary[p, K]
//
// y[p, n]
//
// r[n]
//
// mixN[K]
//
// K[1]
//
// p[1]
//
// n[1]
//
void
ySumBar_j(double* mixsumy,  double* mixbary,  const double* y,  const int* r,  const int* mixN,  const int* K,  
          const int* p,     const int* n);


/***** ***************************************************************************************** *****/
/***** NMix::SS_j                                                                                *****/
/***** ***************************************************************************************** *****/
//
// mixSS[LTp, K]
//
// dwork[p]
//
void
SS_j(double* mixSS,   double* dwork,  const double* mixbary,  const double* y,  const int* r,  const int* K,  
     const int* LTp,  const int* p,   const int* n);


/***** ***************************************************************************************** *****/
/***** NMix::ySum_SSm_j                                                                          *****/
/***** ***************************************************************************************** *****/
//
// mixSSm[LT(p), K]:   mixSSm[,j] = sum_{i: r_i=j}(y_i - mu_j)*(y_i - mu_j)'
//
void
ySum_SSm_j(double* mixsumy,  double* mixSSm,  const double* y,  const int* r,  const double* mu,  const int* K,  
           const int* LTp,   const int* p,    const int* n);


/***** ***************************************************************************************** *****/
/***** NMix::prior_derived                                                                       *****/
/***** ***************************************************************************************** *****/
//
//  Calculate variables derived from the parameters of the prior distribution
//  ===================================================================================================
//
//  /***** logK:                log(1), log(2), ..., log(Kmax)                                                               *****/
//  /***** log_lambda:          log(lambda)                                                                                  *****/
//  /***** c_xi:                c[j]*xi[j], j=0, ..., Kmax-1                                                                 *****/
//  /*****                      * initialize it by xi when priormuQ = MUQ_IC                                                 *****/
//  /***** log_c:               log(c[j]), j=0, ..., Kmax-1                                                                  *****/
//  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
//  /***** sqrt_c:              sqrt(c[j]), j=0, ..., Kmax-1                                                                 *****/
//  /*****                      * initialize it by 0 when priormuQ = MUQ_IC                                                  *****/
//  /***** log_Wishart_const:   Logarithm of the constant in the Wishart density which depends only on degrees of freedom    *****/
//  /***** D_Li:                Cholesky decompositions of D[j]^{-1}, j=0, ..., Kmax-1                                       *****/
//  /*****                      * initialize it by unit matrices when priormuQ = MUQ_NC                                      *****/
//  /***** Dinv_xi:             D[j]^{-1} %*% xi[j], j=0, ..., Kmax-1                                                        *****/
//  /*****                      *initialize it by zero vectors when priormuQ = MUQ_NC                                        *****/
//  /***** log_dets_D:          log_dets based on D matrices                                                                 *****/
//  /*****                      * initialize it by zeros when priormuQ = MUQ_NC                                              *****/
//
//  p[1]
//  priorK[1]
//  priormuQ[1]
//  Kmax[1]
//  lambda[1]
//  xi[p*Kmax]
//  c[Kmax]
//  Dinv[]
//  zeta[1]
//
//  logK[Kmax]
//  log_lambda[1]
//  c_xi[p*Kmax]
//  log_c[Kmax]
//  sqrt_c[Kmax]
//  log_Wishart_const[1]
//  D_Li[LT(p)*Kmax]
//  Dinv_xi[p*Kmax]
//  log_dets_D[2*Kmax]
//
void
prior_derived(const int* p,      
              const int* priorK,  
              const int* priormuQ,  
              const int* Kmax,     
              const double* lambda,  
              const double* xi,  
              const double* c,    
              const double* Dinv,   
              const double* zeta,
              double* logK,  
              double* log_lambda,
              double* c_xi,  
              double* log_c,       
              double* sqrt_c,      
              double* log_Wishart_const,
              double* D_Li,  
              double* Dinv_xi,     
              double* log_dets_D,  
              int*    err);


/***** ***************************************************************************************** *****/
/***** NMix::init_derived                                                                        *****/
/***** ***************************************************************************************** *****/
//
//  Calculate variables derived from the initial values of mixture parameters
//  ===================================================================================================
//
//  /***** log_dets:  log_dets for mixture covariance matrices                                    *****/
//  /***** logw:  Log-weights                                                                     *****/
//  /***** Q:   Mixture inverse variances - compute them from Li                                  *****/
//  /***** Sigma:   Mixture variances - compute them from Li                                      *****/
//  /***** Mean, MeanData:  Mixture overall means                                                 *****/
//  /***** Var, VarData:    Mixture overall variance                                              *****/
//  /***** Corr, CorrData:  Mixture overall std. deviations and correlations                      *****/
//  /***** XiInv:              Diagonal matrix with gamma^{-1}'s on a diagonal                    *****/
//  /***** log_sqrt_detXiInv:  log|XiInv|^{1/2}                                                   *****/    
//
//
// p[1]
// nxw[1]
// Kmax[1]
// K[1]
// distribution[1]     assumed mixture distribution (NORMAL or MVT)
// w[K, nxw]
// mu[p*K]
// Li[LT(p)*K]
// df[K]               initial degrees of freedom for each mixture component (not used when distribution == NMix::NORMAL)
// shift[p]
// scale[p]
// gammaInv[p]
//
// log_dets[Kmax]
// logw[Kmax, nxw]
// Q[LT(p)*Kmax]
// Sigma[LT(p)*Kmax]
// Mean[p, nxw]
// Var[LT(p), nxw]
// Corr[LT(p), nxw]
// MeanData[p, nxw]
// VarData[LT(p), nxw]
// CorrData[LT(p), nxw]
// XiInv[LT(p)]
// log_sqrt_detXiInv[1]
// err[1]
//
void
init_derived(const int* p,    
	     const int* nxw,     
             const int* Kmax,      
             const int* K,  
             const int* distribution,
             const double* w,      
             const double* mu,     
             const double* Li,
             const double* df,
             const double* shift,  
             const double* scale,  
             const double* gammaInv,   
             double* log_dets,  
             double* logw,               
             double* Q,         
             double* Sigma,
             double* Mean,      
             double* Var,                
             double* Corr,
             double* MeanData,  
             double* VarData,            
             double* CorrData,
             double* XiInv,     
             double* log_sqrt_detXiInv,  
             int*    err);


}    /*** end of namespace NMix ***/

#endif
