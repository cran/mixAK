//
//  PURPOSE:   Implementation of methods declared in Dist_TNorm.h
//
// 
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   13/11/2007
//
// ======================================================================
//
#include "Dist_TNorm.h"

namespace Dist{

/***** ***************************************************************************************** *****/
/***** Dist::rTNorm1                                                                             *****/
/***** ***************************************************************************************** *****/
void
rTNorm1(double* x,  const double* mu,  const double* sigma,  const double* a,  const double* b,  const int* trunc)
{
  static double U, Za, Zb, Phia, Phib;

  switch (*trunc){
  case 0:
    Za = (*a - *mu)/(*sigma);
    Phia = pnorm(Za, 0, 1, 1, 0);
    U = Phia + (1 - Phia)*unif_rand();                                    /** U ~ Unif(Phia, 1)                            **/
    if (U > Dist::N_prob1) *x = *mu + *sigma * Dist::N_limit;             /** lower limit is already numerically = Infty   **/
    else{
      if (U < Dist::N_prob0) *x = *mu - *sigma * Dist::N_limit;           /** X = mu + sigma * Phi^{-1}(0)                 **/
      else                   *x = *mu + *sigma * qnorm(U, 0, 1, 1, 0);    /** X = mu + sigma * Phi^{-1}(U)                 **/     
    }
    return;

  case 1:
    *x = *a;
    return;

  case 2:
    Za = (*a - *mu)/(*sigma);
    Phia = pnorm(Za, 0, 1, 1, 0);
    U = Phia*unif_rand();                                                 /** U ~ Unif(0, Phia)                           **/
    if (U < Dist::N_prob0) *x = *mu - *sigma * Dist::N_limit;             /** upper limit is already numerically = -Infty **/
    else{
      if (U > Dist::N_prob1) *x = *mu + *sigma * Dist::N_limit;            /** X = mu + sigma * Phi^{-1}(1)                **/
      else                   *x = *mu + *sigma * qnorm(U, 0, 1, 1, 0);     /** X = mu + sigma * Phi^{-1}(U)                **/     
    }
    return;

  case 3:
    Za = (*a - *mu)/(*sigma);
    Zb = (*b - *mu)/(*sigma);
    Phia = pnorm(Za, 0, 1, 1, 0);         
    Phib = pnorm(Zb, 0, 1, 1, 0);         
    U = Phia + (Phib - Phia)*unif_rand();                                 /** U ~ Unif(Phia, Phib)            **/
    if (U < Dist::N_prob0) *x = *mu - *sigma * Dist::N_limit;             /** X = mu + sigma * Phi^{-1}(0)    **/
    else{
      if (U > Dist::N_prob1) *x = *mu + *sigma * Dist::N_limit;            /** X = mu + sigma * Phi^{-1}(1)    **/
      else                   *x = *mu + *sigma * qnorm(U, 0, 1, 1, 0);     /** X = mu + sigma * Phi^{-1}(U)    **/     
    }
    return;

  case 4:
    *x = *mu + *sigma * norm_rand();
    return;

  default:
    error("Dist::rTnorm1:  Unimplemented value of trunc.\n");
  }
}


#ifdef __cplusplus
extern "C" {
#endif

/***** ***************************************************************************************** *****/
/***** Dist::rTNorm1_R                                                                           *****/
/***** ***************************************************************************************** *****/
void
rTNorm1_R(double* x,      const double* mu,            const double* sigma,          const double* a,  const double* b,  const int* trunc,  
          const int* nx,  const int* mu_sigma_common,  const int* a_b_trunc_common)
{
  int i;
  double *xP;
  const double *muP, *sigmaP, *aP, *bP;
  const int *truncP;

  GetRNGstate(); 
  xP = x;
  if (*mu_sigma_common){
    if (*a_b_trunc_common){
      for (i = 0; i < *nx; i++){
	Dist::rTNorm1(xP, mu, sigma, a, b, trunc);
        xP++;
      }
    }else{
      aP = a;
      bP = b;
      truncP = trunc;
      for (i = 0; i < *nx; i++){
	Dist::rTNorm1(xP, mu, sigma, aP, bP, truncP);
        xP++;
        aP++;
        bP++;
        truncP++;
      }
    }
  }else{
    muP    = mu;
    sigmaP = sigma;
    if (*a_b_trunc_common){
      for (i = 0; i < *nx; i++){
	Dist::rTNorm1(xP, muP, sigmaP, a, b, trunc);
        xP++;
        muP++;
        sigmaP++;
      }
    }else{
      aP = a;
      bP = b;
      truncP = trunc;
      for (i = 0; i < *nx; i++){
	Dist::rTNorm1(xP, muP, sigmaP, aP, bP, truncP);
        xP++;
        aP++;
        bP++;
        truncP++;
        muP++;
        sigmaP++;
      }
    }
  }
  PutRNGstate();

  return;
}

#ifdef __cplusplus
}
#endif


}    /*** end of namespace Dist ***/

