//
//  PURPOSE:   Normal mixture model, birth reversible jump 
//             
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   28/01/2008
//
//  FUNCTIONS:  
//     * RJMCMCbirth   29/01/2008:  Birth of a new component
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_BIRTH_H_
#define _NMIX_RJMCMC_BIRTH_H_

#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>

#include "AK_Basic.h"

#include "Dist_MVN.h"
#include "Dist_Wishart.h"

#include "NMix.h"
#include "NMix_orderComp.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMCbirth                                                                         *****/
/***** ***************************************************************************************** *****/
//
// accept[1]       OUTPUT:  Indicator of acceptance of the proposal
//
// log_AR[1]       OUTPUT:  Logarithm of the acceptance ratio
//
// K[1]            INPUT:   Number of mixture components before the birth move was proposed
//                 OUTPUT:  Number of mixture components after the birth move
//
// w[K+1]          INPUT:   Mixture weights (on the first K places) before the birth move was proposed
//                 OUTPUT:  Mixture weights (on the first K+1 or K places) after the birth move
//
// logw[K+1]       INPUT:   Mixture log-weights (on the first K places) before the birth move was proposed
//                 OUTPUT:  Mixture log-weights (on the first K+1 or K places) after the birth move
//
// mu[p, K+1]      INPUT:   Mixture means (in the first K "columns") before the birth move was proposed
//                 OUTPUT:  Mixture means (in the first K+1 or K "columns") after the birth move
//
// Q[LT(p), K+1]   INPUT:   Mixture inverse variances (in the first K "columns") before the birth move was proposed
//                 OUTPUT:  Mixture inverse variances (in the first K+1 or K "columns") after the birth move
//
// Li[LT(p), K+1]  INPUT:   Cholesky decompositions of the mixture inverse variances (in the first K "columns") before the birth move was proposed
//                 OUTPUT:  Cholesky decompositions of the mixture inverse variances (in the first K+1 or K "columns") after the birth move
//
// Sigma[LT(p), K+1] INPUT:   Mixture variances (in the first K "columns") before the birth move was proposed
//                   OUTPUT:  Mixture variances (in the first K+1 or K "columns") after the birth move
// 
// log_dets[2, K+1]  INPUT:   Factors to evaluate a normal density (in the first K "columns") before the birth move was proposed
//                   OUTPUT:  Factors to evaluate a normal density (in the first K+1 or K  "columns") after the birth move
//
// order[K+1]      INPUT:   Mixture order indeces (on the first K places) before the birth move was proposed
//                OUTPUT:   Mixture order indeces (on the first K+1 or K places) after the birth move 
//
// rank[K+1]        INPUT:   Mixture rank indeces (on the first K plcaes) before the death move was proposed
//                 OUTPUT:   Mixture rank indeces (on the first K+1 or K places) after the birth move 
//
// mixN[K+1]       INPUT:   Numbers of observations belonging to each component
//                 OUTPUT:  Updated numbers of observations belonging to each component (on the first K+1 or K places)
//                          * mixN[j] = sum_{i=1}^n I[r_i=j]
//
// dwork[?]        Working array, see the code for its length
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
// sqrt_c[Kmax]      Square roots of prior precisions of the mixture means when priormuQ = MUQ_NC
//
// log_c[Kmax]       Logarithms of c
// 
// xi[p,Kmax]        Prior means of the mixture means
//
// D_Li[LT(p),Kmax]  Cholesky decompositions of the prior inverse variances for the mixture means when priormuQ = MUQ_IC
//
// log_dets_D[2, Kmax]  log_dets factors related to D matrices
//
// zeta[1]           Prior degrees of freedom of the Wishart distribution on the mixture inverse variances
//
// gammaInv[LT(p)]   Diagonal of the inverted prior scale matrix of the Wishart distribution on the mixture inverse variances
//                   * a priori E(Sigma_j^{-1}) = zeta * Xi and Xi^{-1} = diag(gammaInv)
//
// priormuQ[1]       Type of the prior for the mixture means and (inverse) variances
//                   0 = natural conjugate
//                   1 = independent conjugate
//
void
RJMCMCbirth(int* accept,              double* log_AR,
            int* K,                   double* w,                        double* logw,                double* mu,    
            double* Q,                double* Li,                       double* Sigma,               double* log_dets,  
            int* order,               int* rank,                        int* mixN,
            double* dwork,            int* err,
            const int* p,             const int* n,
            const int* Kmax,          const double* logK,               const double* log_lambda,    const int* priorK,
            const double* logPbirth,  const double* logPdeath,          const double* delta,  
            const double* sqrt_c,     const double* log_c,              const double* xi,            const double* D_Li,               const double* log_dets_D, 
            const double* zeta,       const double* gammaInv,
            const int* priormuQ);

}    /*** end of namespace NMix **/

#endif

