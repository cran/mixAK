//
//  PURPOSE:   Implementation of methods declared in Dist_mixNorm.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   06/11/2008
//
// ======================================================================
//
#include "Dist_mixNorm.h"

namespace Dist{

/***** ********************************************************************************* *****/
/***** Dist::dmixNorm                                                                    *****/
/***** ********************************************************************************* *****/
void
dmixNorm(double* dens,
         const double* x,  
         const int* K,  const double* w,  const double* mu,  const double* sigma)
{
  static int k;
  static const double *wP;
  static const double *muP;
  static const double *sigmaP;

  *dens = 0.0;
  wP = w;
  muP = mu;
  sigmaP = sigma;
  for (k = 0; k < *K; k++){
    *dens += *wP * dnorm(*x, *muP, *sigmaP, 0);
    wP++;
    muP++;
    sigmaP++;    
  }

  return;
}


/***** ********************************************************************************* *****/
/***** Dist::rmixNorm                                                                    *****/
/***** ********************************************************************************* *****/
void
rmixNorm(double* x,   double* dens,
         const int* K,  const double* w,  const double* cumw,  const double* mu,  const double* sigma)
{
  static int r, k;
  static const double *wP;
  static const double *muP;
  static const double *sigmaP;

  /*** Sample the component ***/
  Dist::rDiscrete2(&r, cumw, K);

  /*** Sample the value from the component ***/
  *x = rnorm(mu[r], sigma[r]);

  /*** Evaluate the density ***/
  *dens = 0.0;
  wP = w;
  muP = mu;
  sigmaP = sigma;
  for (k = 0; k < *K; k++){
    *dens += *wP * dnorm(*x, *muP, *sigmaP, 0);
   
    wP++;
    muP++;
    sigmaP++;    
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
dmixNorm_R(double* dens,
           const double* x,   const int* K,     const double* w,  const double* mu,  const double* sigma,
           const int* npoints)
{
  int k, j;
  double *densP;
  const double *xP;

  /*** Evaluate density ***/
  densP = dens;
  xP = x;
  for (j = 0; j < *npoints; j++){
    Dist::dmixNorm(densP, xP, K, w, mu, sigma);    
    densP++;
    xP++;
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
rmixNorm_R(double* x,          double* dens,     double* cumw,
           const int* K,       const double* w,  const double* mu,  const double* sigma,  
           const int* npoints)
{
  int k, j;
  double *cumwP, *xP, *densP;
  const double *wP;

  /*** Computation of cumw                              ***/
  wP    = w;
  cumwP = cumw;

  *cumwP = *wP;
  wP++;
  cumwP++;
  for (k = 1; k < *K; k++){  
    *cumwP = *(cumwP - 1) + (*wP);
    wP++;
    cumwP++;
  }

  /*** Sample and evaluate the density ***/
  GetRNGstate(); 

  densP = dens;
  xP = x;
  for (j = 0; j < *npoints; j++){
    Dist::rmixNorm(xP, densP, K, w, cumw, mu, sigma);

    densP++;
    xP++;
  }

  PutRNGstate();   

  return;
}

#ifdef __cplusplus
}
#endif

}  /** end of namespace Dist **/
