//
//  PURPOSE:   Normal mixture model, split reversible jump 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   03/01/2008
//
//  FUNCTIONS:  
//     * RJMCMCsplit:  Split component into a new one
//              16/01/2008 more or less fully implemented for p <= 2
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_SPLIT_H_
#define _NMIX_RJMCMC_SPLIT_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"
#include "AK_Utils.h"
#include "AK_LAPACK.h"

#include "Dist_MVN.h"
#include "Dist_Wishart.h"

#include "Rand_RotationMatrix.h"

#include "NMix.h"
#include "NMix_orderComp.h"
#include "NMix_RJMCMC_logJac_part3.h"
//#include "NMix_logJacLambdaVSigma.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCsplit: Split move                                                             *****/
/***** ***************************************************************************************** *****/
//
// accept[1]       OUTPUT:  Indicator of acceptance of the proposal
//
// log_AR[1]       OUTPUT:  Logarithm of the acceptance ratio
//
// K[1]            INPUT:   Number of mixture components before the split move was proposed
//                 OUTPUT:  Number of mixture components after the split move
//
// w[K+1]          INPUT:   Mixture weights (on the first K places) before the split move was proposed
//                 OUTPUT:  Mixture weights (on the first K+1 or K places) after the split move
//
// logw[K+1]       INPUT:   Mixture log-weights (on the first K places) before the split move was proposed
//                 OUTPUT:  Mixture log-weights (on the first K+1 or K places) after the split move
//
// mu[p, K+1]      INPUT:   Mixture means (in the first K "columns") before the split move was proposed
//                 OUTPUT:  Mixture means (in the first K+1 or K "columns") after the split move
//
// Q[LT(p), K+1]   INPUT:   Mixture inverse variances (in the first K "columns") before the split move was proposed
//                 OUTPUT:  Mixture inverse variances (in the first K+1 or K "columns") after the split move
//
// Li[LT(p), K+1]  INPUT:   Cholesky decompositions of the mixture inverse variances (in the first K "columns") before the split move was proposed
//                 OUTPUT:  Cholesky decompositions of the mixture inverse variances (in the first K+1 or K "columns") after the split move
//
// Sigma[LT(p), K+1] INPUT:   Mixture variances (in the first K "columns") before the split move was proposed
//                   OUTPUT:  Mixture variances (in the first K+1 or K "columns") after the split move
// 
// log_dets[2, K+1]  INPUT:   Factors to evaluate a normal density (in the first K "columns") before the split move was proposed
//                   OUTPUT:  Factors to evaluate a normal density (in the first K+1 or K  "columns") after the split move
//
// order[K+1]      INPUT:   Order indeces of the mixture components (on the first K places) before the split move was proposed
//                 OUTPUT:  Order indeces of the mixture components (on the first K+1 or K places) after the split move
//
// rank[K+1]       INPUT:   Rank indeces of the mixture components (on the first K places) before the split move was proposed
//                 OUTPUT:  Rank indeces of the mixture components (on the first K+1 or K places) after the split move
//
// r[n]            INPUT:   Allocations before the split move was proposed
//                 OUTPUT:  Updated allocations after the split move
//
// mixN[K+1]       INPUT:   Numbers of observations belonging to each component
//                 OUTPUT:  Updated numbers of observations belonging to each component (on the first K+1 or K places)
//                          * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// rInv[KMax][n]   INPUT:   Indeces of columns of y indicating observations belonging to each component 
//                 OUTPUT:  Updated indeces
//
// u[1 + p + p]    OUTPUT:  Generated values of the auxiliary vector needed for the split move
//                          ASSUMPTION: u1 = u[0] \in (0, 1)
//                                      u2 = u[1], ..., u[p] satisfies: u[p] \in (0, 1) (corresponds to the LARGEST eigen value), u[1], ..., u[p-1] \in (-1, 1)
//                                      u3 = u[p+1], ..., u[2*p] satisfies: u[p+1], ..., u[2*p] \in (0, 1) and u[2*p] corresponds to the LARGEST eigen value
//
// P[p, p]         OUTPUT:   Generated rotation matrix needed for the split move
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
// delta[1]               Prior hyperparameter for the Dirichlet prior on the weights
//
// c[Kmax]                Prior precisions of the mixture means when priormuQ = MUQ_NC
//
// log_c[Kmax]            Logarithms of c
// 
// xi[p,Kmax]             Prior means of the mixture means
//
// D_Li[LT(p),Kmax]       Cholesky decompositions of the prior inverse variances for the mixture means when priormuQ = MUQ_IC
//
// log_dets_D[2, Kmax]    log_dets factors related to D matrices
//
// zeta[1]                Prior degrees of freedom of the Wishart distribution on the mixture inverse variances
//
// log_Wishart_const[1]   Logarithm of the term in the Wishart density which only depends on degrees of freedom
//
// gammaInv[LT(p)]        Diagonal of the inverted prior scale matrix of the Wishart distribution on the mixture inverse variances
//                        * a priori E(Sigma_j^{-1}) = zeta * Xi and Xi^{-1} = diag(gammaInv)
//
// log_sqrt_detXiInv[1]   log|XiInv|^{1/2}
//
// priormuQ[1]            Type of the prior for the mixture means and (inverse) variances
//                        0 = natural conjugate
//                        1 = independent conjugate
//
// pars_dens_u[2*(1 + p + p)]   Parameters of the proposal density of the auxiliary vector u
//
// r_u[FUNCTION]                Function to generate auxiliary vector u and compute its log-density
//
void
RJMCMCsplit(int* accept,           double* log_AR,
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
            void (*r_u)(double* u,  double* log_dens_u,  const double* pars_dens_u,  const int* p));

}

#endif
