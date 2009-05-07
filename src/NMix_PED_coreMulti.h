//
//  PURPOSE:   Normal mixture model, computation of the penalized expected deviance
//             * core part in the multivariate case
//               (part which is the same for uncensored and censored data)
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   08/11/2008
//
//  FUNCTIONS:  
//     * NMix_PED_coreMulti  08/11/2008:  
//
// ====================================================================================================
//
#ifndef _NMIX_PENALIZED_EXPECTED_DEVIANCE_CORE_MULTIVARIATE_H_
#define _NMIX_PENALIZED_EXPECTED_DEVIANCE_CORE_MULTIVARIATE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Basic.h"
#include "Dist_mixMVN.h"

namespace NMix{

/***** ************************************************************************** *****/
/***** NMix::PED_coreMulti                                                        *****/
/***** ************************************************************************** *****/
//
// fy_1[1]              f(y1 | chain1), where y1 is observed value (uncensored observation)
//                      or replicate sampled from the distribution given by the first chain
//                      truncated on the observed interval
//
// fy_2[1]              f(y2 | chain2), where y2 is observed value (uncensored observation)
//                      or replicate sampled from the distribution given by the second chain
//                      truncated on the observed interval
//
// yrep1[p]             replicate from the distribution given by the first chain
//
// yrep2[p]             replicate from the distribution given by the second chain
//
// fyrep1_1[1]          f(yrep1 | chain1)
// 
// fyrep1_2[1]          f(yrep1 | chain2)
// 
// fyrep2_1[1]          f(yrep2 | chain1)
// 
// fyrep2_2[1]          f(yrep2 | chain2)
// 
// pm_indDevObs[1]      += log(f(y1 | chain1)) + log(f(y2 | chain2))
//
// pm_indpopt[1]        += log(f(yrep1 | chain1)/f(yrep1 | chain2)) + log(f(yrep2 | chain2)/f(yrep2 | chain1))
//
// pm_windpopt[1]       += (1/(M*f(y1 | chain1)*f(y2 | chain2))) * pm_indpopt
//
// invalid_indDevObs[1] += 0/1, 1 if f(y1 | chain1) < Dens_ZERO OR f(y2 | chain2) < Dens_ZERO
//                         if 1, then 0 is added to all pm_*
//
// invalid_indpopt[1]   += 0/1, 1 if f(yrep1 | chain1)/f(yrep1 | chain2) < Dens_ZERO OR
//                                   f(yrep2 | chain2)/f(yrep2 | chain1) < Dens_ZERO
//                         if 1, then 0 is added to pm_indpopt and to pm_windpopt
//
// invalid_windpopt[1]  += 0/1, 1 if log(f(y1 | chain1)) + log(f(y2 | chain2)) < EMin
//
// y1[p]                (replicate) of censored observation sampled from the truncated
//                      distribution given by the first chain
//                      or observation itself if it is not censored
// y2[p]                (replicate) of censored observation sampled from the truncated
//                      distribution given by the second chain
//                      or observation itself if it is not censored
//
void
PED_coreMulti(double* fy_1,            double* fy_2,          double* yrep1,          double* yrep2,
              double* fyrep1_1,        double* fyrep1_2,      double* fyrep2_1,       double* fyrep2_2,  
              double* pm_indDevObs,    double* pm_indpopt,    double* pm_windpopt,    double* sum_ISweight,  // double* ch_ISweight,
              int* invalid_indDevObs,  int* invalid_indpopt,  int* invalid_windpopt,  double* work_mix,
              const double* y1,         const int* K1,      const double* w_dets1,  const double* cumw1,
              const double* mu1,        const double* Li1,
              const double* y2,         const int* K2,      const double* w_dets2,  const double* cumw2,  
              const double* mu2,        const double* Li2,
              const int* p,             const int* M,           
              const double *Dens_ZERO,  const double *EMin);

}  /** end of namespace NMix **/

#endif
