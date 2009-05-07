//
//  PURPOSE:   Implementation of methods declared in Dist_MVN.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   05/11/2007
//
// ======================================================================
//
#include "Dist_MVN.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::dMVN1                                                                               *****/
/***** ***************************************************************************************** *****/
void
dMVN1(double* log_dens,  double* work,
      const double* x,   const int* unlog,
      const double* mu,  const double* Li,       const double* log_dets,
      const int* nx,     const int* mu_nonZERO)
{
  static int i;
  static double *dP;
  static const double *cdP1, *cdP2;

  /*** work = x - mu ***/
  if (*mu_nonZERO){
    dP = work;
    cdP1 = x;
    cdP2 = mu;
    for (i = 0; i < *nx; i++){
      *dP = *cdP1 - *cdP2;
      dP++;
      cdP1++;
      cdP2++;
    }
  }else{
    AK_Basic::copyArray(work, x, *nx);
  }

  /*** work = t(Li) %*% (x - mu) ***/
  F77_CALL(dtpmv)("L", "T", "N", nx, Li, work, &AK_Basic::_ONE_INT);          /* Lapack:  work = t(Li) %*% work */

  /*** log_dens = -0.5 * t(x - mu) %*% Li %*% t(Li) %*% (x - mu) ***/
  AK_BLAS::ddot2(log_dens, work, *nx);
  *log_dens *= -0.5;

  /*** log_dens += sum(log(Li[j,j])) - (n/2)*log(2*pi) ***/
  cdP1 = log_dets;
  *log_dens += *cdP1;
  cdP1++;
  *log_dens += *cdP1;
  
  if (*unlog) *log_dens = AK_Basic::exp_AK(*log_dens);
  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::rMVN1                                                                               *****/
/***** ***************************************************************************************** *****/
void
rMVN1(double* x,         double* log_dens,      double* work,
      const double *mu,  const double *Li,      const double *log_dets,           
      const int* nx,     const int* mu_nonZERO)
{
  static int i;
  static double *dP;
  static const double *cdP;

  /*** Sample z ~ N(0, I) ***/
  dP = x;
  for (i = 0; i < *nx; i++){
    *dP = norm_rand();
    dP++;
  }

  /*** Compute -0.5 * t(z) %*% z to get the -0.5 * t(x - mu) %*% Q %*% (x - mu) part of the log-density ***/
  AK_BLAS::ddot2(log_dens, x, *nx);
  *log_dens *= -0.5;

  /*** Solve t(L) %*% v = z,  then v = t(L)^{-1} %*% z ~ N(0, Q^{-1}) ***/
  /*** Store the solution in x                                        ***/
  AK_LAPACK::chol_solve_backward(x, Li, nx);

  /*** Compute x = mu + v, then x = mu + L^{-1}z ~ N(mu, Q^{-1}) ***/
  if (*mu_nonZERO){
    cdP = mu;
    dP = x;
    for (i = 0; i < *nx; i++){
      *dP += *cdP;
      cdP++;
      dP++;
    }
  }

  /*** Add + sum(log(Li[j,j])) - (n/2)*log(2*pi) to the log of the density ***/
  cdP = log_dets;
  *log_dens += *cdP;
  cdP++;
  *log_dens += *cdP;

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::rMVN2                                                                               *****/
/***** ***************************************************************************************** *****/
void
rMVN2(double* x,         double* mu,              double* log_dens,  double* work,
      const double *Li,  const double *log_dets,  const int* nx)
{
  static int i;
  static double *dP;
  static const double *cdP;

  /*** Solve Li %*% w = b, then w = Li^{-1} %*% b  ***/
  /*** Store the solution in mu                    ***/
  AK_LAPACK::chol_solve_forward(mu, Li, nx);

  /*** Solve t(Li) %*% mu = w, then mu = t(Li)^{-1} %*% Li^{-1} %*% b = (Li %*% t(Li))^{-1} %*% b = Q^{-1} %*% b ***/
  AK_LAPACK::chol_solve_backward(mu, Li, nx);


  /***** FROM HERE, IT IS THE SAME AS rMVN1 *****/
  /***** ================================== *****/
  /*** Sample z ~ N(0, I) ***/
  dP = x;
  for (i = 0; i < *nx; i++){
    *dP = norm_rand();
    dP++;
  }

  /*** Compute -0.5 * t(z) %*% z to get the -0.5 * t(x - mu) %*% Q %*% (x - mu) part of the log-density ***/
  AK_BLAS::ddot2(log_dens, x, *nx);
  *log_dens *= -0.5;

  /*** Solve t(L) %*% v = z,  then v = t(L)^{-1} %*% z ~ N(0, Q^{-1}) ***/
  /*** Store the solution in x                                        ***/
  AK_LAPACK::chol_solve_backward(x, Li, nx);

  /*** Compute x = mu + v, then x = mu + L^{-1}z ~ N(mu, Q^{-1}) ***/
  cdP = mu;
  dP = x;
  for (i = 0; i < *nx; i++){
    *dP += *cdP;
    cdP++;
    dP++;
  }

  /*** Add + sum(log(Li[j,j])) - (n/2)*log(2*pi) to the log of the density ***/
  cdP = log_dets;
  *log_dens += *cdP;
  cdP++;
  *log_dens += *cdP;

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::ldMVN1                                                                              *****/
/***** ***************************************************************************************** *****/
void
ldMVN1(double* log_dens,        double* work,
       const double* x,         const double* mu,  const double* Li,       
       const double* log_dets,  const int* nx)
{
  static int i;
  static double *dP;
  static const double *cdP1, *cdP2;

  /*** work = x - mu ***/
  dP = work;
  cdP1 = x;
  cdP2 = mu;
  for (i = 0; i < *nx; i++){
    *dP = *cdP1 - *cdP2;
    dP++;
    cdP1++;
    cdP2++;
  }

  /*** work = t(Li) %*% (x - mu) ***/
  F77_CALL(dtpmv)("L", "T", "N", nx, Li, work, &AK_Basic::_ONE_INT);          /* Lapack:  work = t(Li) %*% work */

  /*** log_dens = -0.5 * t(x - mu) %*% Li %*% t(Li) %*% (x - mu) ***/
  AK_BLAS::ddot2(log_dens, work, *nx);
  *log_dens *= -0.5;

  /*** log_dens += sum(log(Li[j,j])) - (n/2)*log(2*pi) ***/
  cdP1 = log_dets;
  *log_dens += *cdP1;
  cdP1++;
  *log_dens += *cdP1;
  
  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::ldMVN2                                                                              *****/
/***** ***************************************************************************************** *****/
void
ldMVN2(double* log_dens,        double* work,
       const double* x,         const double* mu,  const double* L,       
       const double* log_dets,  const int* nx)
{
  static int i;
  static double *dP;
  static const double *cdP1, *cdP2;

  /*** work = x - mu ***/
  dP = work;
  cdP1 = x;
  cdP2 = mu;
  for (i = 0; i < *nx; i++){
    *dP = *cdP1 - *cdP2;
    dP++;
    cdP1++;
    cdP2++;
  }

  /*** work = L^{-1} %*% (x - mu) ***/
  AK_LAPACK::chol_solve_forward(work, L, nx);

  /*** log_dens = -0.5 * t(x - mu) %*% t(L^{-1}) %*% L^{-1} %*% (x - mu) ***/
  AK_BLAS::ddot2(log_dens, work, *nx);
  *log_dens *= -0.5;

  /*** log_dens += -sum(log(L[j,j])) - (n/2)*log(2*pi) ***/
  cdP1 = log_dets;
  *log_dens += *cdP1;
  cdP1++;
  *log_dens += *cdP1;
  
  return;
}


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::dMVN1_R                                                                           *****/
/***** ***************************************************************************************** *****/
void
dMVN1_R(double* log_dens,  double* Q,              double* work,       int* err,
        const double* x,   const int* unlog,       const double* mu,
        const int* nx,     const int* mu_nonZERO,  const int* npoints)
{
  /*** Cholesky decomposition of Q, i.e., Q = Li %*% t(Li) ***/
  F77_CALL(dpptrf)("L", nx, Q, err);
  if (*err) error("Dist::dMVN1_R: Cholesky decomposition of the precision matrix failed.\n");

  int i;
  double *dP;
  double log_dets[2];
  const double *cdP;

  /*** log_dets[0] = sum(log(Li[j,j])) = log(|Q|^{1/2}) ***/
  cdP = Q;
  dP = log_dets;
  *dP = 0.0;
  for (i = *nx; i > 0; i--){
    *dP += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log_dets[1] = -nx * log(sqrt(2*pi)) ***/
  dP++;
  *dP = -(*nx) * M_LN_SQRT_2PI;

  /*** Evaluate density ***/
  dP = log_dens;
  cdP = x;
  for (i = 0; i < *npoints; i++){
    Dist::dMVN1(dP, work, cdP, unlog, mu, Q, log_dets, nx, mu_nonZERO);
    dP++;
    cdP += *nx;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::rMVN1_R                                                                           *****/
/***** ***************************************************************************************** *****/
void
rMVN1_R(double* x,         double* log_dens,  double* Q,              double* work,        int* err,
        const double *mu,  const int* nx,     const int* mu_nonZERO,  const int* npoints)
{
  /*** Cholesky decomposition of Q, i.e., Q = Li %*% t(Li) ***/
  F77_CALL(dpptrf)("L", nx, Q, err);
  if (*err) error("Dist::rMVN1_R: Cholesky decomposition of the precision matrix failed.\n");

  int i;
  double *dP1, *dP2;
  double log_dets[2];
  const double *cdP;

  /*** log_dets[0] = sum(log(Li[j,j])) = log(|Q|^{1/2}) ***/
  cdP = Q;
  dP1 = log_dets;
  *dP1 = 0.0;
  for (i = *nx; i > 0; i--){
    *dP1 += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log_dets[1] = -nx * log(sqrt(2*pi)) ***/
  dP1++;
  *dP1 = -(*nx) * M_LN_SQRT_2PI;

  /*** Generate random numbers and evaluate the log-density ***/
  GetRNGstate(); 
  dP1 = log_dens;
  dP2 = x;
  for (i = 0; i < *npoints; i++){
    Dist::rMVN1(dP2, dP1, work, mu, Q, log_dets, nx, mu_nonZERO);
    dP1++;
    dP2 += *nx;
  }
  PutRNGstate();

  return;
}


/***** ************************************************************************************************* *****/
/***** Dist::rMVN2_R                                                                                     *****/
/***** ************************************************************************************************* *****/
void
rMVN2_R(double* x,      double* mu,         double* log_dens,       
        double* Q,      double* work,       int* err,
        const int* nx,  const int* npoints)
{
  /*** Cholesky decomposition of Q, i.e., Q = Li %*% t(Li) ***/
  F77_CALL(dpptrf)("L", nx, Q, err);
  if (*err) error("Dist::rMVN2_R: Cholesky decomposition of the precision matrix failed.\n");

  /*** Solve Li %*% w = b, then w = Li^{-1} %*% b  ***/
  /*** Store the solution in mu                    ***/
  AK_LAPACK::chol_solve_forward(mu, Q, nx);

  /*** Solve t(Li) %*% mu = w, then mu = t(Li)^{-1} %*% Li^{-1} %*% b = (Li %*% t(Li))^{-1} %*% b = Q^{-1} %*% b ***/
  AK_LAPACK::chol_solve_backward(mu, Q, nx);

  int i;
  double *dP1, *dP2;
  double log_dets[2];
  const double *cdP;

  /*** log_dets[0] = sum(log(Li[j,j])) = log(|Q|^{1/2}) ***/
  cdP = Q;
  dP1 = log_dets;
  *dP1 = 0.0;
  for (i = *nx; i > 0; i--){
    *dP1 += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log_dets[1] = -nx * log(sqrt(2*pi)) ***/
  dP1++;
  *dP1 = -(*nx) * M_LN_SQRT_2PI;

  /*** Generate random numbers and evaluate the log-density ***/
  GetRNGstate(); 
  dP1 = log_dens;
  dP2 = x;
  for (i = 0; i < *npoints; i++){
    Dist::rMVN1(dP2, dP1, work, mu, Q, log_dets, nx, &AK_Basic::_ONE_INT);
    dP1++;
    dP2 += *nx;
  }
  PutRNGstate();

  return;
}


#ifdef __cplusplus
}
#endif


}  /*** end of namespace Dist ***/

