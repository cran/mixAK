//
//  PURPOSE:   Implementation of methods declared in Stat_BLA.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   19/11/2007
//
// ======================================================================
//
#include "Stat_BLA.h"

namespace Stat{

#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Stat::BLA                                                                                 *****/
/***** ***************************************************************************************** *****/
void
BLA(double* beta,      
    double* sigmaR2,      
    double* L,     
    int*    err,  
    const double* mu,  
    const double* Sigma,  
    const int*    p)
{
  static int i, j, p_1;
  static double betaVbeta;
  static double *betaP, *beta0P, *sigmaR2P;
  static const double *muP;

  *err = 0;
  p_1  = *p - 1;

  betaP    = beta;
  sigmaR2P = sigmaR2;
  for (i = 0; i < *p; i++){
    beta0P = betaP;
    betaP++;

    /*** Compute L[i] = Cholesky decomposition of Sigma[-i,-i] ***/
    AK_BLAS::SPjj(L, betaP, sigmaR2P, Sigma, p, &i);
    F77_CALL(dpptrf)("L", &p_1, L, err);
    if (*err) error("Stat::BLA:  Cholesky decomposition of Sigma[-%d,-%d] failed.\n", i, i);

    /*** Compute beta[1:(p-1),i] = Sigma[-i,-i]^{-1} %*% Sigma[i,-i]                                           ***/
    /*** In the middle, compute also betaVbeta = t(beta[1:(p-1),i]) %*% L[i] %*% t(L[i]) %*% beta[1:(p-1),i]   ***/
    AK_LAPACK::chol_solve_forward(betaP, L, &p_1);
    AK_BLAS::ddot2(&betaVbeta, betaP, p_1);
    AK_LAPACK::chol_solve_backward(betaP, L, &p_1);
    
    /*** Compute beta[0,i] = mu[i] - t(beta[1:(p-1), i]) %*% mu[-i] ***/
    /*** Shift betaP at the same time                               ***/
    *beta0P = 0.0;
    muP     = mu;
    for (j = 0; j < i; j++){
      *beta0P -= *betaP * *muP;
      betaP++;
      muP++;
    }
    *beta0P += *muP;
    muP++;
    for (j = i+1; j < *p; j++){
      *beta0P -= *betaP * *muP;
      betaP++;
      muP++;
    }

    /*** Compute residual variance ***/
    *sigmaR2P -= betaVbeta;
    sigmaR2P++;
  }

  return;
}

#ifdef __cplusplus
}
#endif

}  /*** end of namespace Stat ***/


