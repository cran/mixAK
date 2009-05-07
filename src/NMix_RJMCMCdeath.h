//
//  PURPOSE:   Normal mixture model, death reversible jump 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/01/2008
//
//  FUNCTIONS:  
//     * RJMCMCdeath   29/01/2008:  Death of one of empty components
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_DEATH_H_
#define _NMIX_RJMCMC_DEATH_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "AK_Basic.h"

#include "NMix.h"
#include "NMix_orderComp.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCdeath                                                                         *****/
/***** ***************************************************************************************** *****/
//
// accept[1]       OUTPUT:  Indicator of acceptance of the proposal
//
// log_AR[1]       OUTPUT:  Logarithm of the acceptance ratio
//
// K[1]            INPUT:   Number of mixture components before the death move was proposed
//                 OUTPUT:  Number of mixture components after the death move
//
// w[K]            INPUT:   Mixture weights (on the first K places) before the death move was proposed
//                 OUTPUT:  Mixture weights (on the first K-1 or K places) after the death move
//
// logw[K]         INPUT:   Mixture log-weights (on the first K places) before the death move was proposed
//                 OUTPUT:  Mixture log-weights (on the first K-1 or K places) after the death move
//
// mu[p, K]        INPUT:   Mixture means (in the first K "columns") before the death move was proposed
//                 OUTPUT:  Mixture means (in the first K-1 or K "columns") after the death move
//
// Q[LT(p), K]     INPUT:   Mixture inverse variances (in the first K "columns") before the death move was proposed
//                 OUTPUT:  Mixture inverse variances (in the first K-1 or K "columns") after the death move
//
// Li[LT(p), K]    INPUT:   Cholesky decompositions of the mixture inverse variances (in the first K "columns") before the death move was proposed
//                 OUTPUT:  Cholesky decompositions of the mixture inverse variances (in the first K-1 or K "columns") after the death move
//
// Sigma[LT(p), K]   INPUT:   Mixture variances (in the first K "columns") before the death move was proposed
//                   OUTPUT:  Mixture variances (in the first K-1 or K "columns") after the death move
// 
// log_dets[2, K]    INPUT:   Factors to evaluate a normal density (in the first K "columns") before the death move was proposed
//                   OUTPUT:  Factors to evaluate a normal density (in the first K-1 or K  "columns") after the death move
//
// order[K]        INPUT:   Mixture order indeces before the death move was proposed
//                OUTPUT:   Mixture order indeces (on the first K-1 or K places) after the death move 
//
// rank[K]          INPUT:   Mixture rank indeces before the death move was proposed
//                 OUTPUT:   Mixture rank indeces (on the first K-1 or K places) after the death move 
//
// mixN[K]         INPUT:   Numbers of observations belonging to each component
//                 OUTPUT:  Updated numbers of observations belonging to each component (on the first K-1 or K places)
//                          * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// jempty[K]       OUTPUT:  Indeces of empty components
//
// err[1]          OUTPUT:  Error flag
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
// logPbirth[Kmax]   logarithms of probabilities of the birth move given K
//                   * logPbirth[0]      = log(P(birth | K = 1)) = log(1) = 0
//                   * ...
//                   * logPbirth[Kmax-1] = log(P(birth | K = Kmax)) = log(0) = R_NegInf
//
// logPdeath[Kmax]   logarithms of probabilities of the death move given K
//                     * logPdeath[0]      = log(P(death | K = 1)) = log(0) = R_NegInf
//                     * ...
//                     * logPdeath[Kmax-1] = log(P(death | K = Kmax)) = log(1) = 0
//
// delta[1]          Prior hyperparameter for the Dirichlet prior on the weights
//
void
RJMCMCdeath(int* accept,              double* log_AR,
            int* K,                   double* w,                        double* logw,                double* mu,    
            double* Q,                double* Li,                       double* Sigma,               double* log_dets,  
            int* order,               int* rank,                        int* mixN,
            int* jempty,              int* err,
            const int* p,             const int* n,
            const int* Kmax,          const double* logK,               const double* log_lambda,    const int* priorK,
            const double* logPbirth,  const double* logPdeath,          const double* delta);

}    /*** end of namespace NMix **/

#endif

