//
//  PURPOSE:   Implementation of methods declared in GLMM_fitted_Poisson_Log.h
//
//  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
//             arnost.komarek[AT]mff.cuni.cz
//
//  CREATED:   21/10/2009
//
// ======================================================================
//
#include "GLMM_fitted_Poisson_Log.h"

namespace GLMM{

/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Poisson_Log                                                                  *****/
/*****   PROTOTYPE 1                                                                             *****/
/***** ***************************************************************************************** *****/
void
fitted_Poisson_Log(double* fitted,
                   const double* eta_fixed,  const double* eta_random, 
                   const int* nobs)
{
  static int i;
  static double *fittedP;
  static const double *eta_fixedP, *eta_randomP;

  fittedP     = fitted;
  eta_fixedP  = eta_fixed;
  eta_randomP = eta_random;
  for (i = 0; i < *nobs; i++){
    *fittedP = AK_Basic::exp_AK(*eta_fixedP + *eta_randomP);
    fittedP++;
    eta_fixedP++;
    eta_randomP++;
  }

  return;
}


/***** ***************************************************************************************** *****/
/***** GLMM::fitted_Poisson_Log                                                                  *****/
/*****   PROTOTYPE 2                                                                             *****/
/***** ***************************************************************************************** *****/
void
fitted_Poisson_Log(double* fitted,
                   const double* eta_fixed,  const double* eta_random, 
                   const int* I,             const int* n)
{
  static int i, j;
  static double *fittedP;
  static const double *eta_fixedP, *eta_randomP;
  static const int *nP;

  fittedP     = fitted;
  eta_fixedP  = eta_fixed;
  eta_randomP = eta_random;
  nP          = n; 
  for (i = 0; i < *I; i++){
    for (j = 0; j < *nP; j++){
      *fittedP = AK_Basic::exp_AK(*eta_fixedP + *eta_randomP);
      fittedP++;
      eta_fixedP++;
      eta_randomP++;
    }
    nP++;
  }

  return;
}

}
