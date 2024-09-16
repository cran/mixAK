//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMC_logJacLambdaVSigma.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   09/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMC_logJacLambdaVSigma.h"

namespace NMix {

/***** ************************************************************************************************************************* *****/
/***** NMix::RJMCMC_logJacLambdaVSigma:  log-Jacobian of the Sigma -> (Lambda, V) transformation, Sigma = V %*% Lambda %*% t(V)  *****/
/***** ************************************************************************************************************************* *****/
void
RJMCMC_logJacLambdaVSigma(double* logJac,        double* dlambdaV_dSigma,  double* dwork,            int* iwork,  int* err,
                          const double* Lambda,  const double* V,          const double* Sigma,
                          const int* p,          const int* lambda_descend)
{
  static int LTp, k, j, i, l, lambdaSkip, sign;
  static double *dlambda_dSigmaP, *dV_dSigmaP;
  static double *MPpinv, *work_MPpinvSP, *mijTemp, *mij, *Vdescend;
  static double *MPpinvP, *mijTempP, *mijP, *VdescendP;
  static const double *Vhere, *LambdaStart, *VjP, *ViP, *SigmaP, *LambdaP;

  LTp = (*p * (*p + 1))/2;

  MPpinv        = dwork;
  work_MPpinvSP = MPpinv + LTp;
  mijTemp       = work_MPpinvSP;                   /** mijTemp is of length <= *p - 1 and is nowhere used together with work_MPpinvSP **/
  mij           = work_MPpinvSP + (4 + *p) * *p;
  Vdescend      = mij + (*p - 1) * LTp;            /** needed only if lambda_descend = 1 **/
  // next = Vdescend + *p * *p;                 

  if (*lambda_descend){
    /*** Change the order of columns in V ***/
    VdescendP = Vdescend;
    VjP       = V + (*p - 1) * *p;
    for (j = 0; j < *p; j++){
      ViP = VjP;
      for (i = 0; i < *p; i++){
        *VdescendP = *ViP;
        VdescendP++;
        ViP++;
      }
      VjP -= *p;
    }

    /*** Set up V and pointer lambda order for the rest of the function ***/
    Vhere = Vdescend;
    LambdaStart = Lambda + (*p - 1);
    lambdaSkip = -1;
  }
  else{                      /*** lambda's in ascending order ***/
    Vhere = V;
    LambdaStart = Lambda;
    lambdaSkip = 1;
  }

  /*** Compute dlambda/dSigma ***/
  dlambda_dSigmaP = dlambdaV_dSigma;
  VjP             = Vhere;  
  for (k = 0; k < *p; k++){               /*** loop over eigenvalues ***/
    for (j = 0; j < *p; j++){               /*** loop over columns of Sigma ***/
      *dlambda_dSigmaP = *VjP * *VjP;         /*** dlambda[k]/dSigma[j,j] = V[j,k]^2 ***/
      dlambda_dSigmaP++;
      ViP = VjP + 1;
      for (i = j+1; i < *p; i++){             /*** loop over rows of the lower off-diagonal triangle of Sigma ***/
        *dlambda_dSigmaP = 2 * *VjP * *ViP;   /*** dlambda[k]/dSigma[i,j]= 2 * V[j,k] * V[i,k]                ***/
        dlambda_dSigmaP++;
        ViP++;       
      }
      VjP++;
    }
  }

  /*** Compute dlambda/dV ***/
  dV_dSigmaP = dlambda_dSigmaP;
  LambdaP    = LambdaStart;
  VjP        = Vhere;
  for (k = *p-2; k >=0; k--){                /*** loop over the smallest/largest p-1 eigenvalues ***/

    /*** Compute (lambda[p-k-2]*I - Sigma) and store it in MPpinv ***/
    MPpinvP = MPpinv;
    SigmaP  = Sigma;
    for (j = 0; j < *p; j++){
      *MPpinvP = *LambdaP - *SigmaP;
      MPpinvP++;
      SigmaP++;
      for (i = j+1; i < *p; i++){
        *MPpinvP = - *SigmaP;
        MPpinvP++;
        SigmaP++;
      }
    }

    /*** Compute Moore-Penrose pseudoinverse of (lambda[p-k-2]*I - Sigma) ***/
    AK_LAPACK::MPpinvSP(MPpinv, work_MPpinvSP, err, p);
    if (*err){ 
      Rf_warning("NMix::RJMCMC_logJacLambdaVSigma: Moore-Penrose pseudoinverse failed.\n");
      return;
    }

    /*** Compute m factors, m[[i,j]] = MPpinv * H[[i,j]] * V[,p-k-2], for j=0:(p-1), i=j:(p-1), i.e., there are p*(p+1)/2 m values ***/
    /*** H[[i,j]] is a matrix with h[i,j] = h[j,i] = 1 and zeros otherwise                                                         ***/
    /*** In general, m[[i,j]] is a vector of length p, but we only need rows 0:k                                                   ***/ 
    mijP = mij;
    for (j = 0; j < *p; j++){
      AK_BLAS::SPjxScalar(mijP, MPpinv, VjP, p, &j, &k);             /** V[j,p-k-2] * MPpinv[0:k,j] **/
      mijP += (k + 1);
      ViP = VjP + 1;
      for (i = j+1; i < *p; i++){
        AK_BLAS::SPjxScalar(mijTemp, MPpinv, ViP, p, &j, &k);        /** V[i,p-k-2] * MPpinv[0:k,j] **/
        AK_BLAS::SPjxScalar(mijP, MPpinv, VjP, p, &i, &k);           /** V[j,p-k-2] * MPpinv[0:k,i] **/        
        mijTempP = mijTemp;
        for (l = 0; l <= k; l++){                                    /** V[i,p-k-2] * MPpinv[0:k,j] + V[j,p-k-2] * MPpinv[0:k,i] **/
          *mijP += *mijTempP;
          mijP++;
          mijTempP++;
        }
        ViP++;
      }
      VjP++;
    }

    /*** Extract dV/dSigma ***/
    mijP = mij;
    for (l = 0; l <= k; l++){
      mijTempP = mijP;
      for (j = 0; j < *p; j++){
        for (i = j; i < *p; i++){
          *dV_dSigmaP = *mijTempP;
          mijTempP += (k + 1);
          dV_dSigmaP++;
        }
      }
      mijP++;
    }
    LambdaP += lambdaSkip;   /*** +1 or -1 ***/
  }

  /*** Compute determinant and logarithm of its absolute value ***/
  AK_LAPACK::logDetGE(logJac, &sign, dlambdaV_dSigma, iwork, err, &LTp);  
  if (*err){
    Rf_warning("NMix::RJMCMC_logJacLambdaVSigma: AK_LAPACK::logDet failed.\n");
    return;
  }

  return;
}

}    /*** end of namespace NMix ***/

