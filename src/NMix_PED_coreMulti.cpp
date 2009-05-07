//
//  PURPOSE:   Implementation of methods declared in NMix_PED_coreMulti.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   08/11/2008
//
// ====================================================================================================
//
#include "NMix_PED_coreMulti.h"

namespace NMix{

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
              const double *Dens_ZERO,  const double *EMin)
{
  static double ISweight, ratio1, ratio2, Jtheta1_theta2;

  Dist::dmixMVN(fy_1, work_mix, y1, K1, w_dets1, mu1, Li1, p);
  Dist::dmixMVN(fy_2, work_mix, y2, K2, w_dets2, mu2, Li2, p);

  if (*fy_1 < *Dens_ZERO){
    (*invalid_indDevObs)++;
    (*invalid_indpopt)++;
    (*invalid_windpopt)++;
    // *ch_ISweight = 0.0;
    if (*fy_2 >= *Dens_ZERO){
      *pm_indDevObs += log(*fy_2);
    }
    else{
      (*invalid_indDevObs)++;  
    }
  }
  else{
    ISweight = log(*fy_1);
    if (*fy_2 < *Dens_ZERO){
      (*invalid_indDevObs)++;
      (*invalid_indpopt)++;
      (*invalid_windpopt)++;
      // *ch_ISweight = 0.0;
      *pm_indDevObs += ISweight;
    }
    else{
      ISweight += log(*fy_2);
      *pm_indDevObs += ISweight;

      Dist::rmixMVN(yrep1, fyrep1_1, work_mix, K1, w_dets1, cumw1, mu1, Li1, p);
      Dist::dmixMVN(fyrep1_2, work_mix, yrep1, K2, w_dets2, mu2, Li2, p);

      Dist::rmixMVN(yrep2, fyrep2_2, work_mix, K2, w_dets2, cumw2, mu2, Li2, p);
      Dist::dmixMVN(fyrep2_1, work_mix, yrep2, K1, w_dets1, mu1, Li1, p);

      if (*fyrep1_2 < *Dens_ZERO) ratio1 = *fyrep1_1 / *Dens_ZERO;
      else                        ratio1 = *fyrep1_1 / *fyrep1_2;
      if (*fyrep2_1 < *Dens_ZERO) ratio2 = *fyrep2_2 / *Dens_ZERO;
      else                        ratio2 = *fyrep2_2 / *fyrep2_1;

      // if (*fyrep1_1 < *Dens_ZERO || *fyrep1_2 < *Dens_ZERO || *fyrep2_1 < *Dens_ZERO || *fyrep2_2 < *Dens_ZERO){
      if (ratio1 < *Dens_ZERO || ratio2 < *Dens_ZERO){
        (*invalid_indpopt)++;
        (*invalid_windpopt)++;
        // *ch_ISweight = 0.0;
      }
      else{
        // Jtheta1_theta2 = log(*fyrep1_1) - log(*fyrep1_2) + log(*fyrep2_2) - log(*fyrep2_1);
        Jtheta1_theta2 = log(ratio1) + log(ratio2);
        *pm_indpopt += Jtheta1_theta2;

        if (ISweight < *EMin){
          (*invalid_windpopt)++;
          ISweight = exp(-(*EMin))/(*M);
        }
        else{
          ISweight = exp(-ISweight)/(*M);            /** divide by M to prevent too large numbers in sum_ISweight **/
        }
        *pm_windpopt += ISweight * Jtheta1_theta2;
        // *ch_ISweight = ISweight;
        *sum_ISweight += ISweight;
      }
    }
  }

  return;
}

}  /** end of namespace NMix **/
