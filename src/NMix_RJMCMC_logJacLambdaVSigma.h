//
//  PURPOSE:   Normal mixture model, log-Jacobian of the spectral decomposition d(lambda, V)/d(Sigma),
//             i.e., of the transformation Sigma -> (Lambda, V), where Sigma = V %*% Lambda %*% t(V)
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/01/2008
//
//  FUNCTIONS:  
//     *  11/01/2008:   NMix::logJacLambdaVSigma
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_LOG_JAC_LAMBDA_V_SIGMA_H_
#define _NMIX_RJMCMC_LOG_JAC_LAMBDA_V_SIGMA_H_

#include <R.h>
#include <R_ext/Error.h>

#include "AK_BLAS.h"
#include "AK_LAPACK.h"

namespace NMix{

/***** ************************************************************************************************************************* *****/
/***** NMix::RJMCMC_logJacLambdaVSigma:  log-Jacobian of the Sigma -> (Lambda, V) transformation, Sigma = V %*% Lambda %*% t(V)  *****/
/***** ************************************************************************************************************************* *****/
//
// * Taken partially from the Matlab code of I. Papageorgiou (jacobianVLambda2.m, jacobianVLambda3.m, jacobianVLambda5.m)
//   !!!(I am not sure whether it is theoretically correct)!!!
//
// ASSUMPTION: p >= 2
//
// logJac[1]:                      OUTPUT:  Computed log-Jacobian
// 
// dlambdaV_dSigma[LT(p), LT(p)]:  OUTPUT:  LU factorisation (as returned by LAPACK dgetrf)
//                                          of the matrix of derivatives d(lambda,V)/dSigma
//                      
// dwork[LT(p) + (4+p)*p + (p-1)*LT(p) + p*p]:  Working array
//
// iwork[p]:            Working array for AK_LAPACK::logDet
//
// err[1]:              Error flag
//
// Lambda[p]:           Eigenvalues of Sigma (in ASCENDING order)
//
// V[p, p]:             Eigenvectors of Sigma (in columns)
//
// p[1]:                Dimension
//
// lambda_descend[1]:   If <> 0 then d(lambda,V)/d(Sigma) is computed as if lambda's are stored in 
//                      DESCENDING order
//
void
RJMCMC_logJacLambdaVSigma(double* logJac,        double* dlambdaV_dSigma,  double* dwork,            int* iwork,  int* err,
                          const double* Lambda,  const double* V,          const double* Sigma,
                          const int* p,          const int* lambda_descend);
}

#endif
