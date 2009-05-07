//
//  PURPOSE:   Normal mixture model, combine reversible jump 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   02/01/2008
//
//  FUNCTIONS:  
//     * RJMCMCcombine:  Combine two components into a new one
//              24/01/2008 more or less fully implemented for p <= 2
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_COMBINE_H_
#define _NMIX_RJMCMC_COMBINE_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "AK_LAPACK.h"

#include "Dist_MVN.h"
#include "Dist_Wishart.h"

#include "Rand_SamplePair.h"

#include "NMix.h"
#include "NMix_orderComp.h"
#include "NMix_RJMCMC_logJac_part3.h"
//#include "NMix_logJacLambdaVSigma.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCcombine                                                                       *****/
/***** ***************************************************************************************** *****/
//
// accept[1]       OUTPUT:  Indicator of acceptance of the proposal
//
// log_AR[1]       OUTPUT:  Logarithm of the acceptance ratio
//
// K[1]            INPUT:   Number of mixture components before the combine move was proposed
//                 OUTPUT:  Number of mixture components after the combine move
//
// w[K]            INPUT:   Mixture weights (on the first K places) before the combine move was proposed
//                 OUTPUT:  Mixture weights (on the first K-1 or K places) after the combine move
//
// logw[K]         INPUT:   Mixture log-weights (on the first K places) before the combine move was proposed
//                 OUTPUT:  Mixture log-weights (on the first K-1 or K places) after the combine move
//
// mu[p, K]        INPUT:   Mixture means (in the first K "columns") before the combine move was proposed
//                 OUTPUT:  Mixture means (in the first K-1 or K "columns") after the combine move
//
// Q[LT(p), K]     INPUT:   Mixture inverse variances (in the first K "columns") before the combine move was proposed
//                 OUTPUT:  Mixture inverse variances (in the first K-1 or K "columns") after the combine move
//
// Li[LT(p), K]    INPUT:   Cholesky decompositions of the mixture inverse variances (in the first K "columns") before the combine move was proposed
//                 OUTPUT:  Cholesky decompositions of the mixture inverse variances (in the first K-1 or K "columns") after the combine move
//
// Sigma[LT(p), K]   INPUT:   Mixture variances (in the first K "columns") before the combine move was proposed
//                   OUTPUT:  Mixture variances (in the first K-1 or K "columns") after the combine move
// 
// log_dets[2, K]    INPUT:   Factors to evaluate a normal density (in the first K "columns") before the combine move was proposed
//                   OUTPUT:  Factors to evaluate a normal density (in the first K-1 or K  "columns") after the combine move
//
// order[K]        INPUT:   Order indeces of the mixture components (on the first K places) before the combine move was proposed
//                 OUTPUT:  Order indeces of the mixture components (on the first K-1 or K places) after the combine move
// 
// rank[K]         INPUT:   Rank indeces of the mixture components (on the first K places) before the combine move was proposed
//                 OUTPUT:  Rank indeces of the mixture components (on the first K-1 or K places) after the combine move
//
// r[n]            INPUT:   Allocations before the combine move was proposed
//                 OUTPUT:  Updated allocations after the combine move
//
// mixN[K]         INPUT:   Numbers of observations belonging to each component
//                 OUTPUT:  Updated numbers of observations belonging to each component (on the first K-1 or K places)
//                          * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// rInv[KMax][n]   INPUT:   Indeces of columns of y indicating observations belonging to each component 
//                 OUTPUT:  Updated indeces
//
// u[1 + p + p]    OUTPUT:  Computed values of the auxiliary vector needed for the reversible split move
//
// P[p, p]         OUTPUT:  Computed rotation matrix needed for the reversible split move
//                          (not needed for p = 1 when it is equal to 1)
//
// log_dens_u[1 + 1 + p + p]   Logarithm of the joint density and marginal densities of auxiliary variables u (needed for the acceptance ratio)
//
// dwork[?]        Working array, see the code for its length
// 
// iwork[?]        Working array, see the code for its length
//
// y[p, n]         Observations (in columns)
//
// p[1]            Dimension of the response
//
// n[1]            Number of observations
//
// Kmax[1]           Maximal number of mixture components
//
// logK[Kmax]        logK[j] = log(j+1), j=0,...,Kmax-1
// 
// log_lambda[1]     Logarithm of the prior hyperparameter for the truncated Poisson prior on K
//
// priorK[1]         Type of the prior for K
//                   0 = fixed K
//                   1 = uniform on {1,...,Kmax}
//                   2 = Poisson(lambda) truncated at Kmax
//
// logPsplit[Kmax]   logarithms of probabilities of the split move given K
//                   * logPsplit[0]      = log(P(split | K = 1)) = log(1) = 0
//                   * ...
//                   * logPsplit[Kmax-1] = log(P(split | K = Kmax)) = log(0) = R_NegInf
//
// logPcombine[Kmax]   logarithms of probabilities of the combine move given K
//                     * logPcombine[0]      = log(P(combine | K = 1)) = log(0) = R_NegInf
//                     * ...
//                     * logPcombine[Kmax-1] = log(P(combine | K = Kmax)) = log(1) = 0
//
// delta[1]          Prior hyperparameter for the Dirichlet prior on the weights
//
// c[Kmax]           Prior precisions of the mixture means when priormuQ = MUQ_NC
//
// log_c[Kmax]       Logarithms of c
// 
// xi[p,Kmax]        Prior means of the mixture means
//
// D_Li[LT(p),Kmax]  Cholesky decompositions of the prior inverse variances for the mixture means when priormuQ = MUQ_IC
//
// zeta[1]           Prior degrees of freedom of the Wishart distribution on the mixture inverse variances
//
// log_dets_D[2, Kmax]  log_dets factors related to D matrices
//
// log_Wishart_const[1]   Logarithm of the term in the Wishart density which only depends on degrees of freedom
//
// gammaInv[LT(p)]   Diagonal of the inverted prior scale matrix of the Wishart distribution on the mixture inverse variances
//                   * a priori E(Sigma_j^{-1}) = zeta * Xi and Xi^{-1} = diag(gammaInv)
//
// log_sqrt_detXiInv[1]  log|XiInv|^{1/2}
//
// priormuQ[1]       Type of the prior for the mixture means and (inverse) variances
//                   0 = natural conjugate
//                   1 = independent conjugate
//
// pars_dens_u[2*(1 + p + p)]   Parameters of the proposal density of the auxiliary vector u
//
// ld_u[FUNCTION]    Function to compute log-density of the auxiliary vector u
//
void
RJMCMCcombine(int* accept,           double* log_AR,
              int* K,                double* w,             double* logw,           double* mu,    
              double* Q,             double* Li,            double* Sigma,          double* log_dets,  
              int* order,            int* rank,             int* r,                 int* mixN,         int** rInv,
              double* u,             double* P,             double* log_dens_u,                  
              double* dwork,         int* iwork,            int* err,
              const double* y,          const int* p,                     const int* n,
              const int* Kmax,          const double* logK,               const double* log_lambda,  const int* priorK,
              const double* logPsplit,  const double* logPcombine,        const double* delta,  
              const double* c,          const double* log_c,              const double* xi,          const double* D_Li,               const double* log_dets_D,
              const double* zeta,       const double* log_Wishart_const,  const double* gammaInv,    const double* log_sqrt_detXiInv,  
              const int* priormuQ,      const double* pars_dens_u,
              void (*ld_u)(double* log_dens_u,  const double* u,  const double* pars_dens_u,  const int* p));

}

#endif
