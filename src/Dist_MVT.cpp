//
//  PURPOSE:   Implementation of methods declared in Dist_MVT.h
//
// 
//  AUTHOR:    Arnošt Komárek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  LOG:       20111208  created
//             20161011  lgamma -> lgammafn
//             20220419  FCONE added where needed
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


/***** ***************************************************************************************** *****/
/***** Dist::ldMVT1                                                                               *****/
/***** ***************************************************************************************** *****/
void
ldMVT1(double*       log_dens,
       double*       work,
       const double* x,    
       const double* nu,
       const double* mu,
       const double* Li,
       const double *log_dets,
       const int*    nx)
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
  F77_CALL(dtpmv)("L", "T", "N", nx, Li, work, &AK_Basic::_ONE_INT FCONE FCONE FCONE);          /* Lapack:  work = t(Li) %*% work */

  /*** log_dens = 1 +  (t(x - mu) %*% Li %*% t(Li) %*% (x - mu)) / nu ***/
  AK_BLAS::ddot2(log_dens, work, *nx);
  *log_dens /= *nu;
  *log_dens += 1;

  /*** log_dens = (-(nu + nx)/2) * log(1 +  (t(x - mu) %*% Sigma^(-1) %*% (x - mu)) / nu) ***/
  *log_dens = (-(*nu + *nx) / 2) * log(*log_dens);

  /*** log_dens += log|Sigma|^(-1/2) + log(normal. constant) ***/
  cdP1 = log_dets;
  *log_dens += *cdP1;
  cdP1++;
  *log_dens += *cdP1;

  return;
}


/***** ***************************************************************************************** *****/
/***** Dist::deriv_ldMVT_x                                                                       *****/
/***** ***************************************************************************************** *****/
void
deriv_ldMVT_x(double*       gradient,
              double*       Hessian,
              const double* x,    
              const double* nu,
              const double* mu,
              const double* Q,
              const double* Li,
              const int*    nx)
{
  static int i, j;
  static double one_SS, mult_const1, mult_const2;
  static double *dP1, *dP2;
  static const double *cdP1, *cdP2;

  /*** gradient = x - mu ***/
  dP1 = gradient;
  cdP1 = x;
  cdP2 = mu;
  for (i = 0; i < *nx; i++){
    *dP1 = *cdP1 - *cdP2;
    dP1++;
    cdP1++;
    cdP2++;
  }

  /*** gradient = t(Li) %*% (x - mu) ***/
  F77_CALL(dtpmv)("L", "T", "N", nx, Li, gradient, &AK_Basic::_ONE_INT FCONE FCONE FCONE);          /* Lapack:  gradient = t(Li) %*% gradient */

  /*** one_SS = 1 + (t(x - mu) %*% Q %*% (x - mu)) / nu ***/
  AK_BLAS::ddot2(&one_SS, gradient, *nx);
  one_SS /= *nu;
  one_SS += 1.0;

  /*** gradient = Li %*% t(Li) %*% (x - mu) = Q %*% (x - mu) ***/
  F77_CALL(dtpmv)("L", "N", "N", nx, Li, gradient, &AK_Basic::_ONE_INT FCONE FCONE FCONE);          /* Lapack:  gradient = Li %*% gradient */
  
  /*** Hessian = (2 / (nu * one_SS^2)) * Q %*% (x - mu) %*% t(x - mu) %*% Q ***/
  mult_const1 = 2 / (*nu * one_SS * one_SS);

  dP1   = Hessian;
  cdP1 = gradient;  
  for (j = 0; j < *nx; j++){
    cdP2 = cdP1;
    for (i = j; i < *nx; i++){
      *dP1 = mult_const1 * *cdP1 * *cdP2;
      dP1++;
      cdP2++;
    }
    cdP1++;
  }

  /*** Hessian -= (1 / one_SS) * Q  and then *= (nu + nx) / nu ***/
  /*** gradient *= -(nu + nx) / (nu * one_SS) ***/
  mult_const1 = (*nu + *nx) / *nu;
  mult_const2 = -mult_const1 / one_SS;

  dP1  = Hessian;
  dP2  = gradient;
  cdP1 = Q; 
  for (j = 0; j < *nx; j++){
    *dP2 *= mult_const2;
    dP2++;
    for (i = j; i < *nx; i++){
      *dP1 -= (*cdP1 / one_SS);
      *dP1 *= mult_const1;

      dP1++;
      cdP1++;
    }
  }   

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
  F77_CALL(dpptrf)("L", nx, Q, err FCONE);
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
  *dP1 = lgammafn((*nu + *nx)/2) - lgammafn(*nu / 2) - (*nx) * (0.5 * log(*nu) + M_LN_SQRT_PI);

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


/***** ***************************************************************************************** *****/
/***** Dist::dMVT1_R                                                                           *****/
/***** ***************************************************************************************** *****/
void
dMVT1_R(double*       log_dens,  
        double*       Q,              
        double*       work,       
        int*          err,
        const double* x,   
        const int*    unlog,
        const double* nu,
        const double* mu,      
        const int*    nx,     
        const int*    npoints)
{
  /*** Cholesky decomposition of Q, i.e., Q = Li %*% t(Li) ***/
  F77_CALL(dpptrf)("L", nx, Q, err FCONE);
  if (*err) error("Dist::rMVT1_R: Cholesky decomposition of the precision matrix failed.\n");

  int i;
  double *dP;
  double log_dets[2];
  const double *cdP;

  /*** log_dets[0] = sum(log(Li[j,j])) = log(|Q|^{1/2}) ***/
  cdP = Q;
  dP  = log_dets;
  *dP = 0.0;
  for (i = *nx; i > 0; i--){
    *dP += AK_Basic::log_AK(*cdP);
    cdP += i;
  }

  /*** log_dets[1] ***/
  dP++;
  *dP = lgammafn((*nu + *nx)/2) - lgammafn(*nu / 2) - (*nx) * (0.5 * log(*nu) + M_LN_SQRT_PI);

  /*** Evaluate density ***/
  dP = log_dens;
  cdP = x;
  for (i = 0; i < *npoints; i++){
    Dist::ldMVT1(dP, work, cdP, nu, mu, Q, log_dets, nx);
    if (*unlog) *dP = exp(*dP);
    dP++;
    cdP += *nx;
  }

  return;
}

#ifdef __cplusplus
}
#endif

}  // end of namespace Dist

