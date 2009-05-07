//
//  PURPOSE:   Normal mixture model with reversible jumps,
//             generators and log-densities for auxiliary vector u
//             
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   17/01/2008
//
//  FUNCTIONS:
//      * RJMCMC_r_u_DP  17/01/2008
//
//      * RJMCMC_ld_u_DP  17/01/2008
//
// ====================================================================================================
//
#ifndef _NMIX_RJMCMC_AUX_VECTOR_U_H_
#define _NMIX_RJMCMC_AUX_VECTOR_U_H_

#include <R.h>
#include <Rmath.h>

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::RJMCMC_r_u_DP                                                                       *****/
/***** NMix::RJMCMC_ld_u_DP                                                                      *****/
/***** ***************************************************************************************** *****/
//
// Generating u in the same way as in Dellaportas and Papageorgiou (2006)
// * it is assumed that eigenvalues are sorted in ASCENDING order 
//   (the same as in the Matlab code of I. Papageorgiou, reversal to what's written in their paper)
//
//
// u[1 + p + p]       Generated u vector
//
// log_dens_u[1 + 1 + p + p]  log_dens_u[0] = sum(log_dens_u[1:p])
//                            log_dens_u[1,...] = logarithm of the densities of generated u's
//
// par_dens_u[2* (1 + p + p)]   parameters for the densities of u
//
// p[1]               dimension
//
void
RJMCMC_r_u_DP(double* u,  double* log_dens_u,  const double* pars_dens_u,  const int* p);

void
RJMCMC_ld_u_DP(double* log_dens_u,  const double* u,  const double* pars_dens_u,  const int* p);

}

#endif
