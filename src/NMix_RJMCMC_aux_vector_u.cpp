//
//  PURPOSE:   Implementation of methods declared in NMix_RJMCMC_aux_vector_u.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   17/01/2008
//
// ======================================================================
//
#include "NMix_RJMCMC_aux_vector_u.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMC_r_u_DP                                                                       *****/
/***** ***************************************************************************************** *****/
void
RJMCMC_r_u_DP(double* u,  double* log_dens_u,  const double* pars_dens_u,  const int* p) 
{
  static int i;
  static double *uP, *log_dens_uP;
  static const double *pars_dens_uP;

  uP           = u;
  log_dens_uP  = log_dens_u + 1;
  pars_dens_uP = pars_dens_u;

  /*** ----- u1 ~ Beta(a, b) ----- ***/
  /*** =========================== ***/
  *uP          = rbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  = *log_dens_uP;
  uP++;
  log_dens_uP++;
  pars_dens_uP += 2;

  /*** ----- u2_1,...,u2_{p-1} ~ Unif(a, b) ----- ***/
  /*** ========================================== ***/
  for (i = 1; i < *p; i++){
    *uP          = runif(pars_dens_uP[0], pars_dens_uP[1]);  
    *log_dens_uP = -log(pars_dens_uP[1] - pars_dens_uP[0]);
    *log_dens_u  += *log_dens_uP;
    uP++;
    log_dens_uP++;
    pars_dens_uP += 2;
  }

  /*** ----- u2_p ~ Beta(a, b) ----- ***/
  /*** ============================= ***/
  *uP          = rbeta(pars_dens_uP[0], pars_dens_uP[1]);  
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  += *log_dens_uP;
  uP++;
  log_dens_uP++;
  pars_dens_uP += 2;

  /*** ----- u3_1,...,u3_{p-1} ~ Unif(a, b) ----- ***/
  /*** ========================================== ***/
  for (i = 1; i < *p; i++){
    *uP          = runif(pars_dens_uP[0], pars_dens_uP[1]);  
    *log_dens_uP = -log(pars_dens_uP[1] - pars_dens_uP[0]);
    *log_dens_u  += *log_dens_uP;
    uP++;
    log_dens_uP++;
    pars_dens_uP += 2;
  }

  /*** ----- u3_p ~ Beta(a, b) ----- ***/
  /*** ============================= ***/
  *uP          = rbeta(pars_dens_uP[0], pars_dens_uP[1]);  
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  += *log_dens_uP;
  //uP++;
  //log_dens_uP++;
  //pars_dens_uP += 2;

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::RJMCMC_ld_u_DP                                                                      *****/
/***** ***************************************************************************************** *****/
void
RJMCMC_ld_u_DP(double* log_dens_u,  const double* u,  const double* pars_dens_u,  const int* p) 
{
  static int i;
  static double *log_dens_uP;
  static const double *uP, *pars_dens_uP;

  uP           = u;
  log_dens_uP  = log_dens_u + 1;
  pars_dens_uP = pars_dens_u;

  /*** ----- u1 ~ Beta(a, b) ----- ***/
  /*** =========================== ***/
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  = *log_dens_uP;
  uP++;
  log_dens_uP++;
  pars_dens_uP += 2;

  /*** ----- u2_1,...,u2_{p-1} ~ Unif(a, b) ----- ***/
  /*** ========================================== ***/
  for (i = 1; i < *p; i++){
    *log_dens_uP = -log(pars_dens_uP[1] - pars_dens_uP[0]);
    *log_dens_u  += *log_dens_uP;
    uP++;
    log_dens_uP++;
    pars_dens_uP += 2;
  }

  /*** ----- u2_p ~ Beta(a, b) ----- ***/
  /*** ============================= ***/
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  += *log_dens_uP;
  uP++;
  log_dens_uP++;
  pars_dens_uP += 2;

  /*** ----- u3_1,...,u3_{p-1} ~ Unif(a, b) ----- ***/
  /*** ========================================== ***/
  for (i = 1; i < *p; i++){
    *log_dens_uP = -log(pars_dens_uP[1] - pars_dens_uP[0]);
    *log_dens_u  += *log_dens_uP;
    uP++;
    log_dens_uP++;
    pars_dens_uP += 2;
  }

  /*** ----- u3_p ~ Beta(a, b) ----- ***/
  /*** ============================= ***/
  *log_dens_uP = (pars_dens_uP[0]-1)*log(*uP) + (pars_dens_uP[1]-1)*log(1 - *uP) - lbeta(pars_dens_uP[0], pars_dens_uP[1]);
  *log_dens_u  += *log_dens_uP;
  //uP++;
  //log_dens_uP++;
  //pars_dens_uP += 2;

  return;
}

}

