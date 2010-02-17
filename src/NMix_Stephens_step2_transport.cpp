//
//  PURPOSE:   Implementation of methods declared in NMix_Stephens_step2_transport.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   16/02/2010
//
// =============================================================================
//
#include "NMix_Stephens_step2_transport.h"

namespace NMix{

/***** ***************************************************************************************** *****/
/***** NMix::Stephens_step2_transport                                                            *****/
/***** ***************************************************************************************** *****/
void
Stephens_step2_transport(int*    nchange,
                         int*    order,
                         int*    rank,
                         double* lp_costs,
                         double* lp_solution,
                         int*    lp_r_signs,
                         double* lp_r_rhs,
                         int*    lp_c_signs,
                         double* lp_c_rhs,
                         int*    lp_integers,
                         const double* hatPr_y,
                         const double* Pr_y,
                         const int*    M,
                         const int*    n,
                         const int*    K)
{
  static const double *Pr_yP;

  static int *orderP, *rankP;
  static int m;
  static double minLoss[1];
  static int lp_status[1];
  static int    lp_K[1];
  static int    IZERO[1] = {0};
  static double DZERO[1] = {0};  

  *lp_K = *K;

  lp_costs[0] = 0.0;           /*** start of lp_costs as input to lp_transbig ***/ 

  *nchange = 0;
  orderP   = order;
  rankP    = rank;
  Pr_yP    = Pr_y;
  for (m = 0; m < *M; m++){    /*** loop(m) over number of MCMC iterations ***/

    /***** Calculate the costs matrix *****/
    NMix::Stephens_costMatrix(lp_costs + 1, hatPr_y, Pr_yP, n, K);

    /***** Solve the transportation problem *****/
    //lp_transbig(IZERO, lp_K, lp_K, lp_costs, lp_r_signs, lp_r_rhs, lp_c_signs, lp_c_rhs, minLoss, lp_K, lp_integers, lp_solution, 
    //            IZERO, IZERO, DZERO, DZERO, DZERO, DZERO, DZERO, lp_status);
    // This does not work since lpSolve package does not register its C code to be usable by other packages

    /***** Identify new ranks from supplied solution *****/

  }                            /*** end of loop(m) over number of MCMC iterations ***/  

  return;
}

}    // end of namespace NMix

