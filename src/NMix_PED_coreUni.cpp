//
//  PURPOSE:   Implementation of methods declared in NMix_PED_coreUni.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   08/11/2008
//
// ====================================================================================================
//
#include "NMix_PED_coreUni.h"

namespace NMix{

void
PED_coreUni(double* fy_1,            double* fy_2,          double* yrep1,        double* yrep2,
            double* fyrep1_1,        double* fyrep1_2,      double* fyrep2_1,     double* fyrep2_2,  
            double* pm_indDevObs,    double* pm_indpopt,    double* pm_windpopt,  double* sum_ISweight,  // double* ch_ISweight,
            int* invalid_indDevObs,  int* invalid_indpopt,  int* invalid_windpopt,
            const double* y1,   const int* K1,            const double* w1,  const double* cumw1,  
            const double* mu1,  const double* sigma1,
            const double* y2,   const int* K2,            const double* w2,  const double* cumw2,  
            const double* mu2,  const double* sigma2,
            const int* M,       const double *Dens_ZERO,  const double *EMin)
{
  static double ISweight, ratio1, ratio2, Jtheta1_theta2;

  Dist::dmixNorm(fy_1, y1, K1, w1, mu1, sigma1);
  Dist::dmixNorm(fy_2, y2, K2, w2, mu2, sigma2);

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
 
      Dist::rmixNorm(yrep1, fyrep1_1, K1, w1, cumw1, mu1, sigma1);
      Dist::dmixNorm(fyrep1_2, yrep1, K2, w2, mu2, sigma2);

      Dist::rmixNorm(yrep2, fyrep2_2, K2, w2, cumw2, mu2, sigma2);
      Dist::dmixNorm(fyrep2_1, yrep2, K1, w1, mu1, sigma1);

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
