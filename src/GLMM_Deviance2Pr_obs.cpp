//
//  PURPOSE:   Implementation of methods declared in GLMM_Deviance2Pr_obs.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   29/11/2010
//
// ======================================================================
//
#include "GLMM_Deviance2Pr_obs.h"

extern int iteration;
extern int iter_show;
extern int clus_show;

namespace GLMM{   // start of namespace GLMM

/***** ***************************************************************************************** *****/
/***** GLMM_Deviance2Pr_obs                                                                      *****/
/***** ***************************************************************************************** *****/
void
Deviance2Pr_obs(double* Pr_obs,
                const double* marg_L_i,
                const double* marg_L_ik,
                const double* w,
                const int*    I,               
                const int*    K)
{
  static int i, k;
  static double *Pr_obsP;
  static const double *marg_L_iP, *marg_L_ikP, *wP;

  static double sum_w_L_ik;

  Pr_obsP = Pr_obs;
  marg_L_iP  = marg_L_i;
  marg_L_ikP = marg_L_ik;
  for (i = 0; i < *I; i++){    

    /*** ----- CODE TO CHECK WHETHER marg_L_iP is correct ----- ***/
    //sum_w_L_ik = 0.0;
    //wP = w;
    //for (k = 0; k < *K; k++){
    //  sum_w_L_ik += *wP * *marg_L_ikP;
    //  marg_L_ikP++;
    //  wP++;
    //}
    //marg_L_ikP -= *K;
    /*** ----- END OF CODE TO CHECK WHETHER marg_L_iP is correct ----- ***/    

    wP = w;
    for (k = 0; k < *K; k++){
      *Pr_obsP = (*wP * *marg_L_ikP) / *marg_L_iP;
        
      Pr_obsP++;
      marg_L_ikP++;
      wP++;
    }

    marg_L_iP++;
  }

  return;
}

}

