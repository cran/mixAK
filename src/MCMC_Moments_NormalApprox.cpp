//
//  PURPOSE:   Implementation of methods declared in MCMC_Moments_NormalApprox.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   25/10/2009
//             19/04/2022  FCONE added where needed
//
// ======================================================================
//
#include "MCMC_Moments_NormalApprox.h"

namespace MCMC{

/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox  (PROTOTYPE 1)                                                 *****/
/***** ***************************************************************************************** *****/
void
Moments_NormalApprox(double* cmean, 
                     double* Q,
                     double* log_sqrtdet_Q,
                     double* dwork,
                     int* err,
                     const double* theta,
                     const double* Pprior,
                     const double* P_Mprior,
                     const int* dim,
                     const char* caller)
{
  static int j;
  static double *cmeanP, *QP;
  static const double *PpriorP, *P_MpriorP, *dworkP;

  /*** Canonical mean of the proposal distribution = U + I*theta + Pbeta_Mbeta ***/
  /*** ======================================================================= ***/

    /** Calculate I*theta and store it in dwork **/      
  F77_CALL(dspmv)("L", dim, &AK_Basic::_ONE_DOUBLE, Q, theta, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, dwork, &AK_Basic::_ONE_INT FCONE);

    /** Sum-up U, I*theta, P_Mprior **/
  cmeanP    = cmean;
  dworkP    = dwork;
  P_MpriorP = P_Mprior;
  for (j = 0; j < *dim; j++){
    *cmeanP += *dworkP + *P_MpriorP;
    cmeanP++;
    dworkP++;
    P_MpriorP++;
  }


  /*** Precision matrix of the proposal distribution = I + diag(Pprior) ***/
  /*** ================================================================ ***/
  QP      = Q;
  PpriorP = Pprior;
  for (j = *dim; j > 0; j--){                 /** loop over a diagonal of Li **/
    *QP += *PpriorP;
    QP += j;
    PpriorP++;
  }


  /*** Cholesky decomposition of precision matrix of proposal distribution for theta ***/
  /*** ============================================================================= ***/
  F77_CALL(dpptrf)("L", dim, Q, err FCONE);                 /** this should never fail... **/
  if (*err) Rf_error("%s: Cholesky decomposition of the precision matrix of the proposal distribution failed.\n", caller);


  /*** Compute log(|Q|^{1/2}) = sum(log(Li[j,j])) ***/
  /*** ========================================== ***/
  QP = Q;
  *log_sqrtdet_Q = 0.0;
  for (j = *dim; j > 0; j--){                 /** loop over a diagonal of Li **/
    *log_sqrtdet_Q += AK_Basic::log_AK(*QP);
    QP += j;
  }


  return;
}


/***** ***************************************************************************************** *****/
/***** MCMC::Moments_NormalApprox  (PROTOTYPE 2)                                                 *****/
/***** ***************************************************************************************** *****/
void
Moments_NormalApprox(double* cmean, 
                     double* dwork,                       
                     const double* theta,
                     const double* Imat,
                     const double* P_Mprior,
                     const int* dim)
{
  static int j;
  static double *cmeanP, *QP;
  static const double *P_MpriorP, *dworkP;

  /*** Canonical mean of the proposal distribution = U + I*theta + Pbeta_Mbeta ***/
  /*** ======================================================================= ***/

    /** Calculate I*theta and store it in dwork **/      
  F77_CALL(dspmv)("L", dim, &AK_Basic::_ONE_DOUBLE, Imat, theta, &AK_Basic::_ONE_INT, &AK_Basic::_ZERO_DOUBLE, dwork, &AK_Basic::_ONE_INT FCONE);

    /** Sum-up U, I*theta, P_Mprior **/
  cmeanP    = cmean;
  dworkP    = dwork;
  P_MpriorP = P_Mprior;
  for (j = 0; j < *dim; j++){
    *cmeanP += *dworkP + *P_MpriorP;
    cmeanP++;
    dworkP++;
    P_MpriorP++;
  }

  return;
}


}
