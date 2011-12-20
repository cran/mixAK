//
//  PURPOSE:   Implementation of methods declared in Dist_MVT.h
//
// 
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111208  created
//
// ======================================================================
//
#include "Dist_MVT.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rMVT1                                                                               *****/
/***** ***************************************************************************************** *****/
void
rMVT1(double*       x,
      double*       log_dens,
      const double* nu,
      const double* mu,
      const double* Li,
      const double *log_dets,
      const int*    nx)
{
  static int i;
  static double *xP;
  static const double *muP;
  
  static double v;

  /*** Sample u ~ N(0, I) ***/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /*** Compute t(u) %*% u to get (almost) the t(x - mu) %*% Q %*% (x - mu) part of the log-density. ***/
  /*** It must be multiplied by nu / v at the end.                                                  ***/
  AK_BLAS::ddot2(log_dens, x, *nx);

  /*** Solve t(L) %*% z = u,  then z = t(L)^{-1} %*% u ~ N(0, Q^{-1}) ***/
  /*** Store the solution in x                                        ***/
  AK_LAPACK::chol_solve_backward(x, Li, nx);

  /*** Sample v ~ chisq(nu), resp. v ~ gamma(nu/2, 1/2), for gamma distribution, rate=1/2 -> scale=2 ***/
  v = rgamma(*nu/2, 2.0);   

  /*** Multiply the log-density by nu / v --> numerator and then divide by nu --> fraction in the density expression ***/
  *log_dens /= v;

  /*** Calculate the scale factor of the MVT distribution ***/
  v = sqrt(*nu / v);

  /*** Calculate x = mu + L %*% u * sqrt(nu / v) ***/
  muP = mu;
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP *= v;
    *xP += *muP;
    muP++;
    xP++;
  }

  /*** Finalize calculation of log-density ***/
  *log_dens = (-(*nu + *nx)/2) * log(1 + *log_dens) + log_dets[0] + log_dets[1];
  
  return;
}


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rMVT1_R                                                                           *****/
/***** ***************************************************************************************** *****/
void
rMVT1_R(double* x,         
        double* log_dens,  
        double* Q,              
        int*    err,
        const double* nu,
        const double* mu,  
        const int*    nx,
        const int*    npoints)
{
  /*** Cholesky decomposition of Q, i.e., Q = Li %*% t(Li) ***/
  F77_CALL(dpptrf)("L", nx, Q, err);
  if (*err) error("Dist::rMVT1_R: Cholesky decomposition of the precision matrix failed.\n");

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

  /*** log_dets[1] ***/
  dP1++;
  *dP1 = lgamma((*nu + *nx)/2) - lgamma(*nu / 2) - (*nx) * (0.5 * log(*nu) + M_LN_SQRT_PI);

  /*** Generate random numbers and evaluate the log-density ***/
  GetRNGstate(); 
  dP1 = log_dens;
  dP2 = x;
  for (i = 0; i < *npoints; i++){
    Dist::rMVT1(dP2, dP1, nu, mu, Q, log_dets, nx);
    dP1++;
    dP2 += *nx;
  }
  PutRNGstate();

  return;
}

#ifdef __cplusplus
}
#endif


}  // end of namespace Dist

