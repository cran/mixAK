//
//  PURPOSE:   Implementation of methods declared in NMix_update_sum_Ir_and_sum_Pr_y.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   10/02/2010
//
// =============================================================================
//
#include "NMix_update_sum_Ir_and_sum_Pr_y.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::update_sum_Ir_and_sum_Pr_y                                                          *****/
/***** PROTOTYPE 1                                                                               *****/
/***** ***************************************************************************************** *****/
void
update_sum_Ir_and_sum_Pr_y(int*    sum_Ir,
                           double* sum_Pr_y,
                           double* Pr_y,
                           const double* cum_Pr_y,
                           const int* r,
                           const int* rank,
                           const int* K,
                           const int* n)
{
  static int l, j;
  static int *sum_IrP;
  static double *sum_Pr_yP;
  static const int *rP;
  static const double *Pr_yP;

  AK_Utils::cum_Pr2Pr(Pr_y, cum_Pr_y, K, n);

  rP        = r;
  sum_IrP   = sum_Ir;
  Pr_yP     = Pr_y;
  sum_Pr_yP = sum_Pr_y;

  for (l = 0; l < *n; l++){
    sum_IrP[rank[*rP]]++;
    rP++;
    sum_IrP += *K;

    for (j = 0; j < *K; j++){
      sum_Pr_yP[rank[j]] += *Pr_yP;
      Pr_yP++;
    }
    sum_Pr_yP += *K;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** NMix::update_sum_Ir_and_sum_Pr_y                                                          *****/
/***** PROTOTYPE 2                                                                               *****/
/***** ***************************************************************************************** *****/
void
update_sum_Ir_and_sum_Pr_y(int*    sum_Ir,
                           double* sum_Pr_y,
                           const double* Pr_y,
                           const int* r,
                           const int* rank,
                           const int* K,
                           const int* n)
{
  static int l, j;
  static int *sum_IrP;
  static double *sum_Pr_yP;
  static const int *rP;
  static const double *Pr_yP;

  rP        = r;
  sum_IrP   = sum_Ir;
  Pr_yP     = Pr_y;
  sum_Pr_yP = sum_Pr_y;

  for (l = 0; l < *n; l++){
    sum_IrP[rank[*rP]]++;
    rP++;
    sum_IrP += *K;

    for (j = 0; j < *K; j++){
      sum_Pr_yP[rank[j]] += *Pr_yP;
      Pr_yP++;
    }
    sum_Pr_yP += *K;
  }

  return;
}

}    // end of namespace NMix
