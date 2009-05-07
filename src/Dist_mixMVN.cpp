//
//  PURPOSE:   Implementation of methods declared in Dist_mixMVN.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
// ======================================================================
//
#include "Dist_mixMVN.h"

namespace Dist{

/***** ******************************************************************************** *****/
/***** Dist::dmixMVN                                                                    *****/
/***** ******************************************************************************** *****/
void
dmixMVN(double* dens,      double* work,
        const double* x,  
        const int* K,      const double* w_dets,  
        const double* mu,  const double* Li,        const int* nx)
{
  static int k, i, LTnx;
  static double log_densk;
  static double *workP;
  static const double *xP;
  static const double *w_detsP;
  static const double *muP;
  static const double *LiP;

  LTnx = (*nx * (*nx + 1))/2;

  *dens = 0.0;
  w_detsP = w_dets;
  muP = mu;
  LiP = Li;
  for (k = 0; k < *K; k++){

    /*** work = x - mu ***/
    workP = work;
    xP    = x;
    for (i = 0; i < *nx; i++){
      *workP = *xP - *muP;
      workP++;
      xP++;
      muP++;
    }

    /*** work = t(Li) %*% (x - mu) ***/
    F77_CALL(dtpmv)("L", "T", "N", nx, LiP, work, &AK_Basic::_ONE_INT);          /* Lapack:  work = t(Li) %*% work */

    /*** log_densk = -0.5 * t(x - mu) %*% Li %*% t(Li) %*% (x - mu) ***/
    AK_BLAS::ddot2(&log_densk, work, *nx);
    log_densk *= -0.5;

    /*** *dens += w[k] * |Q|^{1/2} * (2*pi)^(-nx/2) * exp(log_densk) ***/
    *dens += *w_detsP * AK_Basic::exp_AK(log_densk);
    w_detsP++;
    LiP += LTnx;
  }

  return;
}


/***** ******************************************************************************** *****/
/***** Dist::rmixMVN                                                                    *****/
/***** ******************************************************************************** *****/
void
rmixMVN(double* x,         double* dens,          double* work,
        const int* K,      const double* w_dets,  const double* cumw,  
        const double* mu,  const double* Li,      const int* nx)
{
  static int r, k, i, LTnx;
  static double log_densr, log_densk;
  static double *workP;
  static double *xP;
  static const double *w_detsP;
  static const double *muP;
  static const double *LiP;

  LTnx = (*nx * (*nx + 1))/2;

  /*** Sample the component ***/
  Dist::rDiscrete2(&r, cumw, K);


  /*** Sample the value from the component ***/
  w_detsP = w_dets + r;
  muP = mu + *nx * r;
  LiP = Li + LTnx * r;

    /*** Sample z ~ N(0, I) ***/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

    /*** Compute -0.5 * t(z) %*% z to get the -0.5 * t(x - mu) %*% Q %*% (x - mu) part of the log-density ***/
  AK_BLAS::ddot2(&log_densr, x, *nx);
  log_densr *= -0.5;

    /*** Solve t(L) %*% v = z,  then v = t(L)^{-1} %*% z ~ N(0, Q^{-1}) ***/
    /*** Store the solution in x                                        ***/
  AK_LAPACK::chol_solve_backward(x, LiP, nx);

    /*** Compute x = mu + v, then x = mu + L^{-1}z ~ N(mu, Q^{-1}) ***/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP += *muP;
    xP++;
    muP++;
  }

    /*** Evaluate the density contribution of the r-th component ***/
  *dens = *w_detsP * AK_Basic::exp_AK(log_densr);


  /*** Evaluate the density ***/
  w_detsP = w_dets;
  muP = mu;
  LiP = Li;

  k = 0;
  while (k < r){

    /*** work = x - mu ***/
    workP = work;
    xP    = x;
    for (i = 0; i < *nx; i++){
      *workP = *xP - *muP;
      workP++;
      xP++;
      muP++;
    }

    /*** work = t(Li) %*% (x - mu) ***/
    F77_CALL(dtpmv)("L", "T", "N", nx, LiP, work, &AK_Basic::_ONE_INT);          /* Lapack:  work = t(Li) %*% work */

    /*** log_densk = -0.5 * t(x - mu) %*% Li %*% t(Li) %*% (x - mu) ***/
    AK_BLAS::ddot2(&log_densk, work, *nx);
    log_densk *= -0.5;

    /*** *dens += w[k] * |Q|^{1/2} * (2*pi)^(-nx/2) * exp(log_densk) ***/
    *dens += *w_detsP * AK_Basic::exp_AK(log_densk);

    w_detsP++;
    LiP += LTnx;

    k++;
  }

  k++;
  w_detsP++;
  muP += *nx;
  LiP += LTnx;

  while (k < *K){

    /*** work = x - mu ***/
    workP = work;
    xP    = x;
    for (i = 0; i < *nx; i++){
      *workP = *xP - *muP;
      workP++;
      xP++;
      muP++;
    }

    /*** work = t(Li) %*% (x - mu) ***/
    F77_CALL(dtpmv)("L", "T", "N", nx, LiP, work, &AK_Basic::_ONE_INT);          /* Lapack:  work = t(Li) %*% work */

    /*** log_densk = -0.5 * t(x - mu) %*% Li %*% t(Li) %*% (x - mu) ***/
    AK_BLAS::ddot2(&log_densk, work, *nx);
    log_densk *= -0.5;

    /*** *dens += w[k] * |Q|^{1/2} * (2*pi)^(-nx/2) * exp(log_densk) ***/
    *dens += *w_detsP * AK_Basic::exp_AK(log_densk);

    w_detsP++;
    LiP += LTnx;

    k++;
  }

  return;
}


/***** ******************************************************************************** *****/
/***** Dist::dmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
dmixMVN_R(double* dens,      double* w_dets,   double* Li, 
          double* work,      int* err,
          const double* x,   const int* K,     const double* mu,  
          const int* nx,     const int* npoints)
{
  int k, j;
  double wd_tmp;
  double *w_detsP, *LiP, *densP;
  const double *xP;

  const int LTnx = (*nx * (*nx + 1))/2;

  /*** Cholesky decompositions of Q[k], i.e., Q[k] = Li[k] %*% t(Li[k]) ***/
  /*** Computation of w_dets                                            ***/
  w_detsP = w_dets;
  LiP = Li;
  for (k = 0; k < *K; k++){  
    F77_CALL(dpptrf)("L", nx, LiP, err);
    if (*err) error("Dist::dmixMVN_R: Cholesky decomposition of one of the precision matrices failed.\n");

    wd_tmp = -(*nx) * M_LN_SQRT_2PI;
    for (j = *nx; j > 0; j--){
      wd_tmp += AK_Basic::log_AK(*LiP);
      LiP += j;
    }
    *w_detsP *= AK_Basic::exp_AK(wd_tmp);
 
    w_detsP++;
  }

  /*** Evaluate density ***/
  densP = dens;
  xP = x;
  for (j = 0; j < *npoints; j++){
    Dist::dmixMVN(densP, work, xP, K, w_dets, mu, Li, nx);
    densP++;
    xP += *nx;
  }

  return;
}

#ifdef __cplusplus
}
#endif


/***** ******************************************************************************** *****/
/***** Dist::rmixMVN_R                                                                  *****/
/***** ******************************************************************************** *****/
#ifdef __cplusplus
extern "C" {
#endif

void
rmixMVN_R(double* x,         double* dens,     double* w_dets,   double* cumw,  double* Li, 
          double* work,      int* err,
          const int* K,      const double* mu,  
          const int* nx,     const int* npoints)
{
  int k, j;
  double wd_tmp;
  double *w_detsP, *cumwP, *LiP, *xP, *densP;

  /*** Cholesky decompositions of Q[k], i.e., Q[k] = Li[k] %*% t(Li[k]) ***/
  /*** Computation of w_dets                                            ***/
  /*** Computation of cumw                                              ***/
  w_detsP = w_dets;
  cumwP   = cumw;
  LiP = Li;
  for (k = 0; k < *K; k++){  
    F77_CALL(dpptrf)("L", nx, LiP, err);
    if (*err) error("Dist::dmixMVN_R: Cholesky decomposition of one of the precision matrices failed.\n");

    wd_tmp = -(*nx) * M_LN_SQRT_2PI;
    for (j = *nx; j > 0; j--){
      wd_tmp += AK_Basic::log_AK(*LiP);
      LiP += j;
    }
   
    if (k == 0) *cumwP = *w_detsP;
    else        *cumwP = *(cumwP - 1) + *w_detsP;
    *w_detsP *= AK_Basic::exp_AK(wd_tmp);
 
    w_detsP++;
    cumwP++;
  }

  /*** Sample and evaluate the density ***/
  GetRNGstate(); 

  densP = dens;
  xP = x;
  for (j = 0; j < *npoints; j++){
    Dist::rmixMVN(xP, densP, work, K, w_dets, cumw, mu, Li, nx);

    densP++;
    xP += *nx;
  }

  PutRNGstate();   

  return;
}

#ifdef __cplusplus
}
#endif

}   /** end of namespace Dist **/
